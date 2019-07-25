# fate_gpu
This is the gpu code for encryption and decryption in Project Fate

Performance Comparasion (unit: k op/s)

| Device                                                | Encryption | Decryption | Encrypted Addition | Encrypted Multiplication |
|:----------------------------------------------------  |:---------- |:---------- |:------------------ |:------------------------ |
| Intel(R) Xeon(R) Silver 4114 CPU @ 2.20GHz (20 cores) | 9.65       | 34.13      | 234.2              | 143.70                   |
| GTX 1080 Ti                                           | 20.17      | 49.95      | 19327.41           | 1853.54                  |
| Titan Xp                                              | 20.81      | 51.27      | 21588.95           | 1821.53
| Tesla V100                                            | 52.5       | 184.8      | 43516.1            | 6541.1                   |

all test cases are tested with 50000 samples
