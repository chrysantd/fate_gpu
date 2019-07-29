# fate_gpu
This is the gpu code for encryption and decryption in Project Fate

Performance Comparasion (unit: k op/s)

| Device                                                | Encryption | Decryption | Encrypted Addition | Encrypted Multiplication |
|:----------------------------------------------------  |:---------- |:---------- |:------------------ |:------------------------ |
| Intel(R) Xeon(R) Silver 4114 CPU @ 2.20GHz (20 cores) | 9.7        | 34.1       | 234.2              | 143.7                    |
| GTX 1080 Ti                                           | 20.2       | 49.9       | 19327.4            | 1853.5                   |
| Tesla V100                                            | 52.5       | 184.8      | 43516.1            | 6541.1                   |

all test cases are tested with 50000 samples


Notice: for gpu code "pailler.cu", plz use function "fixed_window_powm_odd" to replace "sliding_window_powm_odd" when TPI <= 16