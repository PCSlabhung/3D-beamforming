/*
2025/3/4
this code is to test the accuracy of the approximation of polar coordinate delay calculation
2R/c - (X_channel x cos(phi) x sin(theta) + Y_channel x sin(phi)) / c
phi : elevation angle
theta : angle on planar angle
*/
#include "polar_approximation.h"
#include "T_ord.h"
#include "sin.h"
#include "cos_phi.h"
#include "hls_math.h"
#include "sin_cos_value.h"
#include <ap_fixed.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
//data_t total_delay[33][33][306][64];

void top_model(hls::stream<trans_pkt> &bb_data, hls::stream<trans_pkt> &out_data){
    #pragma HLS interface axis register both port = bb_data name = BUS_IN
    #pragma HLS interface axis register both port = out_data   name = BUS_OUT
    #pragma HLS interface ap_ctrl_none       port = return

    bf_data_t_int local_bb_cos_data[Data_dim][channel_num];
    bf_data_t_int local_bb_sin_data[Data_dim][channel_num];
    copy_mem_v1(bb_data, local_bb_sin_data, local_bb_cos_data);
    top_dataflow(local_bb_cos_data, local_bb_sin_data, out_data);
}
void top_dataflow(bf_data_t_int local_bb_cos_data[Data_dim][channel_num], bf_data_t_int local_bb_sin_data[Data_dim][channel_num], hls::stream<trans_pkt> &out_data){
    #pragma HLS inline off
    #pragma HLS DATAFLOW
    hls::stream<cos_sin_t> cos_stream[channel_num];
    hls::stream<cos_sin_t> sin_stream[channel_num];
    hls::stream<index_t> index_stream[channel_num];
    #pragma HLS STREAM variable= cos_stream     
    #pragma HLS STREAM variable= sin_stream     
    #pragma HLS STREAM variable= index_stream 
    #pragma HLS ARRAY_PARTITION variable = cos_stream dim = 1 type = complete
    #pragma HLS ARRAY_PARTITION variable = sin_stream dim = 1 type = complete
    #pragma HLS ARRAY_PARTITION variable = index_stream dim = 1 type = complete
    #pragma HLS ARRAY_PARTITION variable = local_bb_cos_data dim = 2 type = complete
    #pragma HLS ARRAY_PARTITION variable = local_bb_sin_data dim = 2 type = complete
    polar_approximation(cos_stream, sin_stream, index_stream);
    beamformer(local_bb_cos_data, local_bb_sin_data, cos_stream, sin_stream, index_stream, out_data);

}
void copy_mem(trans_pkt axi_element[Data_dim][channel_num / 2], bf_data_t_int local_sin_cache[Data_dim][channel_num], bf_data_t_int local_cos_cache[Data_dim][channel_num]){
    #pragma HLS inline off
    copy_loop: for(int i = 0; i < Data_dim; i++){
        copy_loop2: for(int j = 0; j < channel_num / 2; j++){
            // in_data_t temp = axi_element[i][j];
            input_data_t temp_value = axi_element[i][j].data;

            local_cos_cache[i][j * 2].range()     = temp_value.range(15, 0);
            local_cos_cache[i][j * 2 + 1].range() = temp_value.range(31, 16);
            local_sin_cache[i][j * 2].range()     = temp_value.range(47, 32);
            local_sin_cache[i][j * 2 + 1].range() = temp_value.range(63, 48);
            #ifndef __SYNTHESIS__
                bf_data_t_int temp1 = local_cos_cache[i][j * 2];
                bf_data_t_int temp2;
                temp2.range() = temp_value.range(15, 0);
                if(temp1 != temp2)
                    cout << "before transformation : " << temp1 << "after transformation : " << temp2 << "\n";
            #endif
    
        }
    }
}
void copy_mem_v1(hls::stream<trans_pkt> &axi_element, bf_data_t_int local_sin_cache[Data_dim][channel_num], bf_data_t_int local_cos_cache[Data_dim][channel_num]){
    #pragma HLS inline off
    copy_loop:for(int i = 0; i < Data_dim; i++){
        copy_loop_2:for(int j = 0; j < channel_num / 2 ;j++){
            trans_pkt temp_value = axi_element.read();
            local_cos_cache[i][j * 2].range()     = temp_value.data.range(15, 0);
            local_cos_cache[i][j * 2 + 1].range() = temp_value.data.range(31, 16);
            local_sin_cache[i][j * 2].range()     = temp_value.data.range(47, 32);
            local_sin_cache[i][j * 2 + 1].range() = temp_value.data.range(63, 48);
            #ifndef __SYNTHESIS__
                bf_data_t_int temp1 = local_cos_cache[i][j * 2];
                bf_data_t_int temp2;
                temp2.range() = temp_value.data.range(15, 0);
                if(temp1 != temp2)
                    cout << "before transformation : " << temp1 << "after transformation : " << temp2 << "\n";
            #endif
        }
    }
}
void beamformer(bf_data_t_int bb_cos_data[Data_dim][channel_num], bf_data_t_int bb_sin_data[Data_dim][channel_num], hls::stream<cos_sin_t> cos_stream[channel_num], hls::stream<cos_sin_t> sin_stream[channel_num], hls::stream<index_t> index_stream[channel_num], hls::stream<trans_pkt> &out_data){
    
    bf_data_t rrf_data_sample_real_cos[channel_num];
    bf_data_t rrf_data_sample_real_sin[channel_num];
    bf_data_t rrf_data_sample_real[no_lines_phi][no_lines][NR];
    bf_data_t rrf_data_sample_imag[no_lines_phi][no_lines][NR];
    bf_data_t CF;
    bf_data_t CF_real, CF_imag;
    #pragma HLS ARRAY_PARTITION variable = rrf_data_sample_real_cos     dim = 1 type = complete
    #pragma HLS ARRAY_PARTITION variable = rrf_data_sample_real_sin     dim = 1 type = complete
	#pragma HLS ARRAY_PARTITION variable = bb_cos_data                  dim = 2 type = complete
    #pragma HLS ARRAY_PARTITION variable = bb_sin_data                  dim = 2 type = complete

    // reset_loop: for(int i = 0; i < no_lines_phi; i++){
    //     reset_loop2: for(int j = 0; j < no_lines; j++){
    //         #pragma HLS pipeline
    //         reset_loop3: for(int k = 0; k < NR; k++){
    //             #pragma HLS pipeline
    //             DATA[i][j][k] = 0;
    //             rrf_data_sample_imag[i][j][k] = 0;
    //         }
    //     }
    // }
    #ifndef __SYNTHESIS__
        cout << "begin beamform\n";          
    #endif
    approx_loop: for(int i = 0; i < no_lines_phi; i++){
		approx_loop2: for(int j = 0; j < no_lines; j++){ // phi
            // #ifndef __SYNTHESIS__
            //     cout << "begin beamform " << i << " " << j << "frame\n";          
            // #endif
			approx_loop3: for(int k = 0; k < NR; k++){ // theta
                #pragma HLS pipeline
				approx_loop4: for(int l = 0; l < channel_num; l++){
                    #pragma HLS unroll
                    index_t index_temp = index_stream[l].read();
                    cos_sin_t cos_phase = cos_stream[l].read();
                    cos_sin_t sin_phase = sin_stream[l].read();
                    bf_data_t_int cos_value = bb_cos_data[index_temp][l];
                    bf_data_t_int sin_value = bb_sin_data[index_temp][l];

                    rrf_data_sample_real_cos[l] = cos_value * cos_phase - sin_value * sin_phase;
                    rrf_data_sample_real_sin[l] = sin_value * cos_phase + cos_value * sin_phase;
				}

                // real_part_sum: for(int l = 0; l < channel_num; l++){
                //     #pragma HLS unroll
                //     DATA[j][k][i] += rrf_data_sample_real_cos[l];
                // }
                rrf_data_sample_real[i][j][k] = sum_pipeline(rrf_data_sample_real_cos);
                // imag_part_sum: for(int l = 0; l < channel_num; l++){
                //     #pragma HLS unroll
                //     rrf_data_sample_imag[j][k][i] -= rrf_data_sample_real_sin[l];
                // }
                rrf_data_sample_imag[i][j][k] = sub_pipeline(rrf_data_sample_real_sin);
                CF_real = 0;
                CF_imag = 0;

                // cos_sqr_sin_sqr: for(int l = 0; l < channel_num; l++){
                //     #pragma HLS unroll
                //     CF_real += rrf_data_sample_real_cos[l] * rrf_data_sample_real_cos[l];
                //     CF_imag += rrf_data_sample_real_sin[l] * rrf_data_sample_real_sin[l];
                // }
                CF_real = MAC_pipeline(rrf_data_sample_real_cos, rrf_data_sample_real_cos);
                CF_imag = MAC_pipeline(rrf_data_sample_real_sin, rrf_data_sample_real_sin);
                if(CF_real + CF_imag == 0){
                    CF = 0;
                }
                else{
                    CF = (rrf_data_sample_real[i][j][k] * rrf_data_sample_real[i][j][k] + rrf_data_sample_imag[i][j][k] * rrf_data_sample_imag[i][j][k]) / (CF_real + CF_imag);
                }
                // out_data[i][j][k].real = rrf_data_sample_real[i][j][k] * CF;
                // out_data[i][j][k].imag  = rrf_data_sample_imag[i][j][k] * CF;
                bf_data_t_64 out_data_temp;
                // out_data_temp.real = rrf_data_sample_real[i][j][k] * CF;
                // out_data_temp.imag  = rrf_data_sample_imag[i][j][k] * CF;
                out_data_temp = CF * CF * (rrf_data_sample_real[i][j][k] * rrf_data_sample_real[i][j][k] + rrf_data_sample_imag[i][j][k] * rrf_data_sample_imag[i][j][k]);
                trans_pkt out_temp;
                out_temp.data.range() = out_data_temp.range();
                out_temp.last = (i == no_lines_phi - 1 && j == no_lines - 1 && k == NR - 1) ? 1 : 0;
                out_data.write(out_temp);
			}
            // #ifndef __SYNTHESIS__
            //     cout << "write to stream finish scan line i = " << i << " j = " << j << endl; 
            // #endif
		}
	}
}
bf_data_t sum_pipeline(bf_data_t to_sum_element[channel_num]){
#pragma HLS inline off
    bf_data_t sum[channel_num] = {0};
    sum[0] = to_sum_element[0];
	#pragma HLS array_partition variable=sum dim = 0 type = complete
	#pragma HLS ARRAY_PARTITION variable=to_sum_element dim = 0 type=complete
    #pragma HLS pipeline II = 1
    for(int i = 1; i < channel_num; i++){
		#pragma HLS unroll
        sum[i] = sum[i - 1] + to_sum_element[i];
    }
    return sum[channel_num - 1];
}
bf_data_t sub_pipeline(bf_data_t to_sum_element[channel_num]){
#pragma HLS inline off
bf_data_t sum[channel_num] = {0};
	#pragma HLS array_partition variable=sum dim = 0 type = complete
	#pragma HLS ARRAY_PARTITION variable=to_sum_element dim = 0 type=complete
    
    sum[0] = -to_sum_element[0];
    #pragma HLS pipeline II = 1
    for(int i = 1; i < channel_num; i++){
        #pragma HLS unroll
        sum[i] = sum[i - 1] - to_sum_element[i];
    }
    return sum[channel_num - 1];
}
bf_data_t MAC_pipeline(bf_data_t to_sum_element[channel_num], bf_data_t to_mul_element[channel_num]){
	bf_data_t sum[channel_num] = {0};
    #pragma HLS array_partition variable=sum dim = 0 type = complete
	#pragma HLS ARRAY_PARTITION variable=to_sum_element dim = 0 type=complete
	#pragma HLS ARRAY_PARTITION variable=to_mul_element dim = 0 type=complete
    #pragma HLS pipeline II = 1
    sum[0] = to_sum_element[0] * to_mul_element[0];
    for(int i = 1; i < channel_num; i++){
        #pragma HLS unroll
        sum[i] = sum[i - 1]  + to_mul_element[i] * to_sum_element[i];
    }
    return sum[channel_num - 1];

}
void polar_approximation(hls::stream<cos_sin_t> cos_stream[channel_num], hls::stream<cos_sin_t> sin_stream[channel_num],hls::stream<index_t> index_stream[channel_num]){
	
    
	/*
	R/c - (X_channel * cos(phi) * sin(theta) + y_channel * sin(phi)) / c
	*/

    

	#pragma HLS inline off
    #pragma HLS ARRAY_PARTITION variable = chx                  type = complete
    #pragma HLS ARRAY_PARTITION variable = chy                  type = complete
    #pragma HLS ARRAY_PARTITION variable = cos_phi              type = complete
    #pragma HLS ARRAY_PARTITION variable = sin_pt               type = complete
    #pragma HLS ARRAY_PARTITION variable = T_ord        dim = 1 type = complete

    // for local cache data_flow
    //copy_mem(bb_cos_data, local_bb_cos_data);
    //copy_mem(bb_sin_data, local_bb_sin_data);

   
    #ifndef __SYNTHESIS__
        cout << "begin delay index calculation\n";          
    #endif

	approx_loop: for(int i = 0; i < no_lines_phi; i++){
		approx_loop2: for(int j = 0; j < no_lines; j++){ // phi
			approx_loop3: for(int k = 0; k < NR; k++){ // theta
                #pragma HLS pipeline
				approx_loop4: for(int l = 0; l < channel_num; l++){
                    #pragma HLS unroll
                    data_t T_ord_temp = T_ord[k][l];
                    data_t T_sd_rd = static_cast<data_t>(chx[l]) * cos_phi[i] * sin_pt[j] + static_cast<data_t>(chy[l]) * sin_pt[i];
                    data_t total_delay    = (T_ord_temp - T_sd_rd) * 1000 / 340;
                    data_t phase_temp  	  = total_delay * 80 * PI;
					cos_sin_t cos_phase   = approx_cos(phase_temp);
					cos_sin_t sin_phase   = approx_sin(phase_temp);
                    //cos_sin_t cos_phase = hls::cos(phase_temp);
                    //cos_sin_t sin_phase = hls::sin(phase_temp);
					index_t index_temp    = total_delay * static_cast<data_t>(192);
                    
                    cos_stream[l].write(cos_phase);
                    sin_stream[l].write(sin_phase);
                    index_stream[l].write(index_temp);
                    
				}
                
			}
		}
	}
	/*
	X_scan = R * cos(phi) * sin(theta);
	Y_scan = R * sin(phi);
	Z_scan = R * cos(phi) * cos(theta);
	*/
	
	
	return ;
}


cos_sin_t approx_sin(data_t theta) {
    // 将 theta 归一化到 [0, 2π) ==> [0, 1)
    // theta = fmodf(theta, TWO_PI);
    // if (theta < 0)
    //     theta += TWO_PI;

	#pragma HLS LATENCY  max = 2
	#pragma HLS PIPELINE  II = 1
    data_t t = theta * static_cast<data_t> (0.159154943092);
    t = t - (int)t;
    cos_sin_t result;

    int quadrant = (int)(t * 4); // 得到象限：0,1,2,3
    data_t angle;
    int sign = 1;
    switch (quadrant) {
        case 0:
            angle = t * 4; // 映射到 [0,1]
            break;
        case 1:
            angle = 2 - t * 4; // 映射到 [0,1]
            break;
        case 2:
            angle = t * 4 - 2; // 映射到 [0,1]
            sign = -1;
            break;
        case 3:
            angle = 4 - t * 4; // 映射到 [0,1]
            sign = -1;
            break;
    }
    int index = (int)(angle * (LUT_SIZE - 1));
    result = sign * sin_cos_lut[index];

    return result;
}

// 根据输入角度 theta（单位：弧度），利用 LUT 和对称性查找余弦值
cos_sin_t approx_cos(data_t theta) {
    // 利用余弦的性质： cos(theta) = sin(theta + π/2)
    return approx_sin(theta + static_cast<data_t>(PI / 2));
}


// void polar_approximation(data_t output[33][33][306][64]){
	
// chagpt version

// cos_sin_t approx_sin(data_t theta) {
//     #pragma HLS inline
//     #pragma HLS PIPELINE II=1

//     // theta ∈ [0, 2π)，我們預設它已經被映射好（不含負值）
//     // normalize 到 [0,1)：theta / 2π
//     // 若你有 LUT_FULL_SCALE_FACTOR 可直接傳進來，減少除法

//     const data_t scale = static_cast<data_t>(LUT_SIZE / (M_PI * 2));
//     int full_index = (int)(theta * scale * 4);  // 0 ~ 4*LUT_SIZE-1

//     int quadrant = full_index >> 10;            // 2 bit: 0 ~ 3
//     int base_index = full_index & 0x3FF;        // 10 bit index

//     index_t lut_index;
//     bool is_neg = false;

//     switch (quadrant) {
//         case 0:
//             lut_index = base_index;
//             is_neg = false;
//             break;
//         case 1:
//             lut_index = LUT_SIZE - 1 - base_index;
//             is_neg = false;
//             break;
//         case 2:
//             lut_index = base_index;
//             is_neg = true;
//             break;
//         case 3:
//             lut_index = LUT_SIZE - 1 - base_index;
//             is_neg = true;
//             break;
//     }

//     sincos_t val = sin_cos_lut[lut_index];
//     return is_neg ? -val : val;
// }
