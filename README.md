# Polar Approximation

## Design Files
- `polar_approximation.cpp`
- `polar_approximation.h`

## Testbench Files
- `polar_approximation_tb.cpp`
- `bb_real.txt`
- `bb_imag.txt`

## Input
- **Channels**: 64  
- **Samples per channel**: 2000  
- **Format**: Already demodulated baseband data  

## Output
- **Dimensions**: `33 × 33 × 324` scan points  

## Scan Range
- **Horizontal (azimuth)**: `-60° ~ +60°`  
- **Vertical (elevation)**: `-60° ~ +60°`  
- **Range**: `0.3 m ~ 1.5 m`  

## Interface
- **I/O**: AXI Stream  

---

### Notes
- Input data (`bb_real.txt`, `bb_imag.txt`) should be provided in text format.  
- The design implements polar approximation for beamforming.  

