These folders contain all of the finalized code required to interface and operate the LSM6DSO IMU, LIS3MDL Magnetometer, MS5607-02BA03 Altimeter, Digi XBee 3 Pro, and Radio Frequency Front End boards.

This interface is intended to be performed with a Raspberry Pi 5B running Raspbian Lite 64-bit.

I2C and UART peripherals are required for this code to function, as well as particular pinouts for the GPIO specified in the attached schematic.

This repository DOES NOT include code required to interface the ADALM-PLUTO Rev C SDRs with the Raspberry Pi 5B.

Code for interfacing the ADALM-PLUTO Rev C SDRs with the Raspberry Pi 5B can be found at:

https://github.com/SL-UAV-MQP/pi_C_code2 
https://github.com/SL-UAV-MQP/pi_C_code

All code in this repository and related repositories is a portion of materials for a Major Qualifying Project (MQP), submitted to the faculty of Worcester Polytechnic Institute (WPI) in partial fulfillment of the requirements for the Degree of Bachelor of Science in Electrical and Computer Engineering as completed by John Frahm, Carthene McTague, Ann Phan, Heath Sainio, and John Song.
