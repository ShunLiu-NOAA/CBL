# CBL
09/22/2020
1. The QVP is used to estimate convective boundary layer depth using ZDR observation at ~4.5 deg elevation angle.
2. create github repository test test test test
3. replied Dave's email.o

09/18/2020
Hi Shun,

 

Attached is code in fortran 90 for calculating a QVP from a time series of 2d ZDR at 4.5 degrees elevation and then calculating a time series of CBL depth by finding the minimum value in the vertical column for each observation time.  I have added what I hope are a lot of comments.  The code compiles, but I have no way of testing it in the current format. 

 

When you are ready and have time, please take a look and ask questions. 

 

Enjoy the weekend!

09/11/2020
Dave:

Hi Shun,

I hope you are your family are healthy and doing fine.  Now that classes have started, meaning that all the work we’ve done in preparation is beginning to relax a little bit, I actually have some time for research again.  I’ve looked through the radar code you sent, but admit that I haven’t found what I am looking for.  I am hoping that you can point me to the right place in the code.
 

I was thinking of providing a subroutine to calculate the QVP and return it.  Is there a place in the current code where a 2-d array of ZDR is available from the same radar scan with dimensions of radials and azimuth angles?  For the QVP we typically use the ~4.5 degree elevation angle, as it is available in both clear-air and precipitation scanning modes.  I would also need the distance from the radar for all the range gates so that I can calculate the height of the beam above the ground.  Please let me know. 
 
Many thanks, and take good care.

Shun's reply
David,

Glad to hear from you. I don't start to work on CBL code yet since the project I am working on is postonded from the end of August to the end of this month. But merging your CBL package to the current radar data flow is in m mind. I will start working on this as soon as possible once my current project pass onsite acceptance test at the end of this month. 

Yes, 2D ZDR at 4.5 degree  elevation is available from radar data flow. If you can provide a subroutine or a single along package to test your algorithm using ZDR at 4.5 elevation, it will speed up the processing. In this case, you can start your test with ZDR at 4.5 with your test dataset and once the subroutine or package is ready, I will dump ZDR from radar processing package using the data format you used and merge your subroutine or package to radar data processing. From my side, I can put ZDR at 4.5, azimuth angle of each beam and distance from the radar for all range gates to a file for each qualified volume scan. your subroutine or package will process this file and generate CBL.

Again, sorry for any inconvenience due to the delay.


