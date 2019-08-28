# LTE_RACH_positioning
source code for UE localization based on LTE RACH signals

your comments are valuable for me and I will be happy to discuss them. The related paper is in online:
A. Fedorov, H. Zhang and Y. Chen, "User Localization Using Random Access Channel Signals in LTE Networks with Massive MIMO," 2018 27th International Conference on Computer Communication and Networks (ICCCN), Hangzhou, 2018, pp. 1-9.
doi: 10.1109/ICCCN.2018.8487359
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8487359&isnumber=8487301

1. The UE location is set in Scenario_plain.m as pUe and then put into UE.loc structure.
2. The time delay between BS and UE is T_LoS in the LoS path and T_NLoS in a NLoS path. 
________________________________________________________________________________________
3. Since the LTE signal has a sampling rate (i.e. the signal is discrete in the time domain), the main feature of modeling the time delay is that the time delay is split into discrete shifts and non-discrete shifts. A discrete shift indicates how many discrete samples the transmitted signal should be moved, and a non-discrete shift indicates how much phase shift should be applied to the discreetly shifted signal. 

!!!IMPORTANT!!! the simulation of time delay consists of two steps: (a) discrete shift without phase rotation and (b) phase rotation according to the non-discrete shift. 
________________________________________________________________________________________
4. LoS and NLoS signals are add together to get a multipath channel effect.
5. Signal processing. 
