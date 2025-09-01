#include "polar_approximation.h"
#include "sin.h"
#include "cos_phi.h"
#include <ap_fixed.h>
#include <math.h>
#include <stdint.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()

const float pi = 3.14159265358979323846f;
const float c = 340;
data_t delay_approx[33][33][306][64];
data_t phase_hw[33][33][306][64];
// cos_sin_t sin_hw[33][33][306][64];
// cos_sin_t cos_hw[33][33][306][64];
index_t index_hw[33][33][306][64];

bf_data_t_int bb_cos_data[Data_dim][channel_num];
bf_data_t_int bb_sin_data[Data_dim][channel_num];
// bf_data_t bf_img[no_lines_phi][no_lines][NR];
// bf_data_t bf_real[no_lines_phi][no_lines][NR];
out_complex_t hw_out[no_lines_phi][no_lines][NR];
in_data_t input_data[Data_dim][channel_num / 2];

int main(){
    float angle[33];
    float range[306];
    hls::stream<trans_pkt> bb_data; 
    hls::stream<trans_pkt> out_data;
    bool test = 0;
    cout << "Test bench started" << endl;

    cout<< "Memory allocation done" << endl;

    ifstream bb_cos_file, bb_sin_file;
    ofstream bf_real_file_out, bf_imag_file_out, bf_abs_file_out, bf_abs_file_out_2;
    bb_cos_file.open("/home/hung52852/beamforming/3D/top/bb_real.txt");
    bb_sin_file.open("/home/hung52852/beamforming/3D/top/bb_imag.txt");
  
    if(bb_cos_file.is_open()){
        cout << "bb_cos_data.txt opened" << endl;
    }
    else{
        cout << "bb_cos_data.txt not opened" << endl;
    }
    if(bb_sin_file.is_open()){
        cout << "bb_sin_data.txt opened" << endl;
    }
    else{
        cout << "bb_sin_data.txt not opened" << endl;
    }
    for(int i = 0; i < Data_dim; i++){
        for(int j = 0; j < channel_num; j++){
           bb_cos_file >> bb_cos_data[i][j];
           bb_sin_file >> bb_sin_data[i][j];
          
            if(bb_cos_data[i][j] != 0 || bb_sin_data[i][j] != 0){
                cout << "bb_cos_data[" << i << "][" << j << "] = " << bb_cos_data[i][j] << endl;
                cout << "bb_sin_data[" << i << "][" << j << "] = " << bb_sin_data[i][j] << endl;
            }
        }
    }
    bb_cos_file.close();
    bb_sin_file.close();
    cout << "bb_cos_data.txt and bb_sin_data.txt read done" << endl;
    ofstream sin1, sin0, cos1, cos0, input_hex;
    // sin1.open("sin1.txt");
    // sin0.open("sin0.txt");
    // cos1.open("cos1.txt");
    // cos0.open("cos0.txt");
    // input_hex.open("input_hex.txt");
    for(int i = 0; i < Data_dim; i++){
        for(int j = 0; j < channel_num / 2; j++){
            input_data_t in;
            in.range(15, 0)  = bb_cos_data[i][j * 2].range();
            in.range(31, 16) = bb_cos_data[i][j * 2 + 1].range();
            in.range(47, 32) = bb_sin_data[i][j * 2].range();
            in.range(63, 48) = bb_sin_data[i][j * 2 + 1].range();
            trans_pkt in_stream;          
            // in_stream.data.range(15, 0)  = bb_cos_data[i][j * 2].range();
            // in_stream.data.range(31, 16) = bb_cos_data[i][j * 2 + 1].range();
            // in_stream.data.range(47, 32) = bb_sin_data[i][j * 2].range();
            // in_stream.data.range(63, 48) = bb_sin_data[i][j * 2 + 1].range();
            if(in != 0){
                //cout << "in type input success" << endl;
            }
            in_stream.data.range() = in.range();
            ap_uint<64> check, check2;
            check.range() = in_stream.data.range();
            check2.range() = in.range();
            // if(check != check2){
            //     cout << "fail input" << endl;
            //     return 1;
            // }
            // if(check != 0){
            //     cout << "no zero" << endl;
            // }
            // if(i == Data_dim - 1 && j == channel_num/2 - 1){
            //     in_stream.last = 1;
            // }
            // else{
            //     in_stream.last = 0;
            // }
            bb_data.write(in_stream);
            // cout << "tb write to stream" << in_stream.data << endl;
        }
    }
    top_model(bb_data, out_data);
    cout << "Approximation done" << endl;
    sin1.close();
    sin0.close();
    cos1.close();
    cos0.close();
    input_hex.close();


    // float angle_diff = (2 * pi / 3) / 32;
    // for(int i = 0; i < 33; i++){
    //     angle[i] = -pi / 3 + i * angle_diff;
    // }


    // ////////////////////////////////////////////////////////////////////
    // // 這邊的range原本是0.2~1.5，但fieldII的matlab file是0.3~1.5
    // // 所以這邊range的值改成是0.3~1.5
    // ////////////////////////////////////////////////////////////////////
    // float z_diff = (1.5 - 0.3) / 305;
    // for(int i = 0; i < 306; i++){
    //     range[i] = 0.3 + i * z_diff;
    // }

    ////////////////////////////////////////////////////////////////////
    // 用cmath裡的sin cos直接去算會跟matlab答案不一樣，
    // 因為matlab裡的sin_th = -1*(no_lines-1)/2*d_th  + (jj-1)*d_th;
    // 不是算angle去查表
    ////////////////////////////////////////////////////////////////////

    // bf_real_file_out.open("bf_real.txt");
    // bf_imag_file_out.open("bf_imag.txt");
    // if(bf_real_file_out.is_open()){
    //     cout << "bf_real.txt opened" << endl;
    // }
    // else{
    //     cout << "bf_real.txt not opened" << endl;
    // }
    // if(bf_imag_file_out.is_open()){
    //     cout << "bf_imag.txt opened" << endl;
    // }
    // else{
    //     cout << "bf_imag.txt not opened" << endl;
    // }
    bf_abs_file_out.open("bf_abs.txt");
    bf_abs_file_out_2.open("bf_abs_hex.txt");
    int count = 0;
    true_value_loop:
    for(int i = 0; i < no_lines_phi; i++){
        for(int j = 0; j < no_lines; j++){ // phi
            for(int k = 0; k < NR; k++){ // theta
                // bf_real_file_out << hw_out[i][j][k].real << endl;
                // bf_imag_file_out << hw_out[i][j][k].imag << endl;
                bf_data_t_64 to_file;
                trans_pkt out_stream;
                out_stream = out_data.read();
                to_file.range() = out_stream.data.range();
                int to_hex = static_cast<int>(to_file.range());
                bf_abs_file_out << to_file << endl;
                bf_abs_file_out_2 << std::hex << std::uppercase<< std::setfill('0') << std::setw(16) << out_stream.data.range() << endl;
                cout << "i = " << i << "j = " << j << "k = "<< k << endl;
            }
        }
        cout << "finish "<< i << " frame \n";
    }
    cout << "True value calculation done" << endl;
    bf_real_file_out.close();
    bf_imag_file_out.close();

    return 0;
}
