###################################
#USRPN210 - use this for an additional EMI detector 
#
#Author: Sai Abitha Srinivas 
#Modified last : 3/10/2023
#
#Libraries needed to run - UHD, Basic Python libraries 
#
##################################

####################################
# Add all imports and headers here
####################################
import uhd
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import time

##########################
# Set all variables here
##########################
num_samps 	= 2560 		# number of samples received
center_freq 	= 2.07338e6 	# Hz
sample_rate 	= 0.2e6 	# Hz
gain 		= 6 		# dB
bw 		= 20e3		# Hz
sleep_time	= 0.1 		# Sleep time to allow locking on refclk
num_phase 	= 67		# Number of phase encodes
num_avg		= 1		# Number of averages
num_par		= 1		# Number of slices/partitions
num_samp_blocks = num_phase * num_avg * num_par

###########################
# Initialize and set USRP
###########################
usrp = uhd.usrp.MultiUSRP()

usrp.set_rx_rate(sample_rate, 0)
usrp.set_rx_freq(uhd.libpyuhd.types.tune_request(center_freq), 0)
usrp.set_rx_gain(gain, 0)
usrp.set_clock_source("external")
time.sleep(sleep_time)		# Wait to allow locking of external refclk
usrp.set_time_source("_external_")
#usrp.set_rx_bandwidth(bw)

###################################################
# Initialize and setup stream and receiver buffer
###################################################
st_args = uhd.usrp.StreamArgs("fc64", "sc16")
st_args.channels = [0]
metadata = uhd.types.RXMetadata()
streamer = usrp.get_rx_stream(st_args)
#recv_buffer = np.zeros((1, num_samps*2), dtype=np.complex128)

#########################################################
# Get Rx gate signal from PPS input and start streaming 
#########################################################
samples = np.zeros(num_samp_blocks*num_samps, dtype=np.complex128)
for i in range(num_samp_blocks):
    recv_buffer = np.zeros((1, num_samps), dtype=np.complex128)
    #Get Rx gate signal from PPS
    time_at_last_pps = usrp.get_time_last_pps().get_real_secs()
    while time_at_last_pps == usrp.get_time_last_pps().get_real_secs():
        time.sleep(0.1) 
    #print("Clock locked?   ", usrp.get_mboard_sensor("ref_locked", 0))
    #print(time_at_last_pps)
    #Start streaming
    stream_cmd = uhd.types.StreamCMD(uhd.types.StreamMode.start_cont)
    stream_cmd.stream_now = True
    streamer.issue_stream_cmd(stream_cmd)
    time.sleep(0.05) 
    # Receive Samples
    streamer.recv(recv_buffer, metadata)
    samples[i*num_samps:(i+1)*num_samps] = recv_buffer[0]
    print("Sample count: ",i)
    print(samples[i*num_samps:(i+1)*num_samps])
    time.sleep(0.2)
    # Stop Stream
    stream_cmd = uhd.types.StreamCMD(uhd.types.StreamMode.stop_cont)
    streamer.issue_stream_cmd(stream_cmd)

########################
# Write out .mat file
########################
time_stamp=timestr = time.strftime("%m%d_%H%M%S")
path="/Users/abithasrinivas/Dropbox/Demo_lowfield_studygroup/3rdchannel_data/"
outFile="{}{}.mat" .format(path,time_stamp)
output={"sample":samples}
savemat(outFile,output)




