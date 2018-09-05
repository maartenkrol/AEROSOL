      program spuv_v2
*****6******************************************************************
* This program applies Langley-Analysis on SPUV-data. 
* It creates a new  
* dataset which does not contain unuseful points which correspond to
* cloud passages.
*
* Original version: Otto Hasekamp, April 1998
*
* Improvements by Bas Henzing in 2000.
*
* Continuous improvement by Piet Stammes, 1998-2003.
*
*
* Modifications:
* 25 August 1998: option to read files from a file list.
* 4 June 1999: Extended to include processing of 1999 and 2000 data.
*              This version can now process data from 1997 onwards.
* 30 July 1999: inclusion of sigma (standard dev.) of *b_aer* in the 
*               calculation of alfa.
* 24 Jan. 2000: Changing the airmass formula.
* 16 Feb. 2001: Extended to include processing of 2001 data.
* 20 Feb. 2001: Corrected 2000 leap year error in *clcdate*.
* 01 Mar. 2002: Extended to include processing of 2002 data.
* 18 Mar. 2002: Bug in routine *angstrom* corrected.
* 10 Oct. 2002: Values of ozone and Rayleigh optical thickness improved.
* 30 Dec. 2002: Modifications related to printing in case of errors.
*               This version yields output in files *channelX.out* 
*               which is ready to be 
*               quality-screened by programme *screen.f*.
* 10 Jan. 2003: Extended to include processing of 2002 data.
* 13 Jan. 2003: Print yearfraction yyyy.xxx in outputfile channelX.dat.
* 18 Feb. 2003: Bug corrected in the cloud filter of subroutine 
*               firstfilter, by setting nsum to 1.
* 25 Feb. 2003: Read and use the auxiliary data file readaux.out, which
*               contains the actual pressure and ozone column data.
*
*****6******************************************************************
*
* Explanation of arrays and variables:
*-------------------------------------
* 
* The arrays below belong to one day of SPUV measurements. 
* The first dimension, with range 1-2 (array(2,....)), indicates
* the morning (1) or afternoon (2) part of the data.
* The variable m=1,2 is used for this AM/PM distinction.
*
* nn(2,nchan): total number of data points per half day per channel.
* nl(2,nchan): number of useful data points per half day per channel,
*                after applying the filter criteria.
*
* The decision whether a day can still be used for Langley
* analysis is set by the threshold *proc_limit*. This is the minimum
* percentage of useful points needed for Langley analysis.
*
* The decision whether a day can still be used for deriving the Angstrom
* exponent is set by the threshold *nuse_limit*. This is the minimum
* number of wavelengths for which b_aer must be available from 
* the Langley analysis.
*****6******************************************************************

	implicit none

	integer i,j,m,n,maxsamples,dimmax,nmax,nsum,nfil,r,daynr,totfil
	integer yr,nmonth,nyear,nday,nhour,nmin,nsec
	integer nn(2),alln(2),nuse(2),nchan,index,nb_aer(2),norg(2)
        integer nuse_limit 
	parameter(nuse_limit=3) 
	parameter(nchan=6)
	parameter(maxsamples=1500)
	parameter(dimmax=5000)
	integer nl(2,nchan)
	integer ierr_lang(2,nchan)
	integer crit(2,dimmax,nchan),ntot(2),wrkchan(nchan)
	real x(2,dimmax,nchan),y(2,dimmax,nchan),xlon,xlat
	real xorg(2,dimmax,nchan),yorg(2,dimmax,nchan)
	real err_oz(nchan), toterr(2,nchan)
	real od(2,dimmax,nchan),fluxfac
        real year,yearfrac
	real allx(2,dimmax),ally(2,dimmax,nchan)
	real proc(2,nchan),dayfr(2,dimmax),lambda(nchan),alfa(2),chi2(2)
	real a(2,nchan),siga(2,nchan),b(2,nchan),sigalfa(2)
	real sigb(2,nchan),sd(2,nchan)
	real b_aer(2,nchan),I_ex(2,nchan)
	real defbral(nchan),defb_ozon(nchan)
	real bral(nchan),b_ozon(nchan)
	real SPC,SPD,SPA(7),SPB(7)
        real proc_limit 
	parameter(proc_limit=5.0) 
	integer odef,pdef 
	real defozcol,defpsurf,ozcol,psurf
	character*10 fname
	character*29 path1
	character*2 path2
	character*2 yrchar
	character*4 yeardir
	character*20 allfname
	character*20 outfname
	logical ex
	parameter(xlon=-5.18)
	parameter(xlat=52.10)

	defpsurf=1013.00
	defozcol=334.0

*****6*****************************************************************
* The path of the dayfiles must be specified here.
* The years 1997-2002 are on the /nobackup disk.
* The current year is on the backup disk.
*****6******************************************************************
        path1='/nobackup/users/stammes/SPUV/'
        path2='./'

*****6*****************************************************************
* Below the wavelengths of the SPUV channels are given.
*****6******************************************************************

	data lambda /368.5, 501.4, 674.6, 780.4, 871.2, 938.9/

*****6******************************************************************
* In array *wrkchan* the working channels are specified: 
* wrkchan(ichannel) = 0/1: not working/working.
*
* When channels are added or removed at a certain date this array 
* has to be changed. See the start of the "do 1000" loop.
* 
* The measurements in Channels 2 and 5 before 7 July 1997 were 
* saturated, so they were not useful.
* These channels were removed for repair until 20 Jan. 1998.
* So, before 20 Jan 1998 channels 2 and 5 were in effect absent.
* After 20 Jan 1998 all 6 channels were again available.
*
*****6******************************************************************
* Below the values for ozone absorption and Rayleigh scattering
* optical thickness, as calculated by the DAK program for a mid 
* latitude
* summer, are set. These values are scaled according
* to measurements of the ozone column and ground pressure
*
*****6******************************************************************
* Note: the default (i.e. standard atmosphere) values of 
* b_ozon and bral were updated on 10 Oct. 2002, to be in agreement
* with the latest DAK calculations for an MLS atmosphere.
* (the previous values of b_ozone (as mentioned in Hasekamp's report
* on p. 32, Table 2), were differing by about 0.002 at 500 nm 
* and 0.001 at 870 nm.
*****6******************************************************************
      data defb_ozon /0.00008,0.01191,0.01388,0.00255,0.00000,0.00000/

      data defbral /0.50868,0.14205,0.04242,0.02352,0.01508,0.01106/

*****6******************************************************************
*  Read the input file spuv.in.
*
*  If index=1 the filenames are listed in a file called ' filelist ',
*  created by the Unix command ' ls -1 '.
*  In that case the ground pressure and ozone column are assumed to be 
*  default (i.e. MLS standard atmosphere values).
*  If index=0 the filenames and ozone and pressure values are read from
*  the file ' spuv.in '.
*
*****6******************************************************************
      open (unit=5, file='spuv.in')
      open (unit=51, file='filelist')

      do i=1,6
	   read(5,'(a)')
      enddo

      read (5,*) index
      read (5,*)       
      read (5,*) totfil

      write(*,900) index, totfil
     
  900 format(i2,'=index: 0=read lines; 1= read filenames from filelist',
     & /,i5,'=number of SPUV dayfiles that have to be analyzed ')

      do 1000 nfil=1,totfil
         do i=1,nchan
            wrkchan(i)=1
         enddo

         if (index.eq.1) then
            read (51,'(a)') fname
         else 
	    read (5,'(a)') fname
         endif
	    write(*,'(3a)') 'Reading file ',fname

         yrchar = fname(1:2)

         if (yrchar.eq.'97') then
*****6******************************************************************
* Before 20 Jan 1998 channels 2 and 5 were in effect absent.
*****6******************************************************************
            wrkchan(2)=0
            wrkchan(5)=0
            yr=97
            year=1997.0
	    inquire(file=path1//'DAYFILES/1997/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/1997/'//fname)
            endif
         endif

         if (yrchar.eq.'98') then
            yr=98
            year=1998.0
	    inquire(file=path1//'DAYFILES/1998/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/1998/'//fname)
            endif
         endif

         if (yrchar.eq.'99') then
            yr=99
            year=1999.0
	    inquire(file=path1//'DAYFILES/1999/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/1999/'//fname)
            endif
         endif

         if (yrchar.eq.'00') then
            yr=00
            year=2000.0
	    inquire(file=path1//'DAYFILES/2000/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/2000/'//fname)
            endif
         endif

         if (yrchar.eq.'01') then
            yr=01
            year=2001.0
	    inquire(file=path1//'DAYFILES/2001/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/2001/'//fname)
            endif
         endif

         if (yrchar.eq.'02') then
            yr=02
            year=2002.0
	    inquire(file=path1//'DAYFILES/2002/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path1//'DAYFILES/2002/'//fname)
            endif
         endif

         if (yrchar.eq.'03') then
            yr=03
            year=2003.0
	    inquire(file=path2//'DAYFILES/2003/'//fname,exist=ex)
	    if (.not.ex) then
	       write(*,'(3a)') 'File ',fname,' does not exist!'
               stop
            else
	       open(1,file=path2//'DAYFILES/2003/'//fname)
            endif
         endif
*****6******************************************************************
* If index=1, the filenames of the dayfiles to be analysed are read from
* file *filelist*. The corresponding ozone column and surface pressure
* are read from file *readaux.out*. 
*
* If index=0, the filenames of the dayfiles to be analysed are read from
* the input file *spuv.in*, as well as the ozone column and pressure.
*
*  default ozone column:     YES: odef=1,  NO: odef=0
*  default surface pressure: YES: pdef=1,  NO: pdef=0
*
*****6******************************************************************

         if (index.eq.1) then
            call auxdata(fname,defozcol,defpsurf,ozcol,psurf)
         else 
            read(5,*) odef
            read(5,*) pdef
            if (odef.eq.0) then
               read(5,*) ozcol
            else
               ozcol=defozcol
            endif
            if (pdef.eq.0) then
               read(5,*) psurf
            else
               psurf=defpsurf
            endif
         endif

*****6******************************************************************
* Scale the ozone and Rayleigh optical thicknesses according to the 
* chosen values for ozcol and psurf (see Eqs. 43-44, report O. Hasekamp)
*
* The values of defbral and defb_ozone are the defaults, which should 
* not be altered! (this was a bug discovered on 25 Feb. 2003, P.S.).
*
*****6******************************************************************
	 do i=1,nchan
	    b_ozon(i)=defb_ozon(i)*(ozcol/defozcol)
            err_oz(i)=0.04*b_ozon(i)
	    bral(i)=defbral(i)*(psurf/defpsurf)
c            write(*,*)'b_ozon',b_ozon(i)
c            write(*,*)'bral',bral(i)
	 enddo

*****6******************************************************************
* Read the data.
*****6******************************************************************
	   call readf(x,y,nn,xorg,yorg,norg,daynr,dayfr,nmax,maxsamples,
     &     ntot,alln,allx,SPC,SPD,SPA,SPB,fname)
     
	   call clcdate(daynr,yr,nmonth,nday,year,yearfrac)

	   call wrdayvariation(y,x,dayfr,nmonth,nday,nn,nl,
     &     od,fluxfac,fname,m,i)

*****6******************************************************************
* In each subroutine where unuseful datapoints are rejected, only the
* remaining data points are stored.
*****6******************************************************************
	   do m=1,2
              nuse(m)=0
	      do i=1,6
	         if(wrkchan(i).eq.1)then

	            call firstcheck(x,y,nl,m,i,fname)

		    call firstfilter(x,y,nmax,nl,nsum,m,i)
 
		    call secondfilter(x,y,nl,nsum,nmonth,nday,m,i)

		    call thirdfilter(x,y,nl,m,i)

	            if (ntot(m).gt.0) then
	               proc(m,i)=(float(nl(m,i))/float(ntot(m)))*100.
	            else
	               proc(m,i)=0
	            endif

*****6******************************************************************
* Langley analysis is only performed if the 
* percentage of useful points is at least *proc_limit*.
*
* The index *ierr_lang(m,i)* indicates if a Langley regression was
* successful for that channel:
* ierr_lang=0: successful; no error 
* ierr_lang=1: negative b_aer; error on btot is larger, so acceptable
* ierr_lang=1: negative b_aer; error on btot is smaller, so unacceptable
* ierr_lang=99: no Langley regression was performed 
*
* Find the number of useful Langley regression channels, *nuse*, 
* to determine the Angstrom relation for this day 
* (the water vapour channel, channel 6, cannot be used).
*****6******************************************************************
		    if (proc(m,i).ge.proc_limit) then
		       call langley(x,y,nl,fname,daynr,a,
     &                      siga,b,sigb,b_aer,sd,I_ex,fluxfac,
     &                      b_ozon,err_oz,bral,m,i,toterr,ierr_lang)
                    else
                       ierr_lang(m,i)=99
                    endif

                    if ((ierr_lang(m,i).le.1).and.(i.ne.6)) then 
                       nuse(m)=nuse(m)+1
		    endif
	         endif
	      enddo

*****6******************************************************************
* The Angstrom parameter is only derived if the number of channels for 
* which *b_aer* is derived is at least *nuse_limit*.
*****6******************************************************************
c DEBUG:
c              write (*,902) m, nuse(m)
c  902 format('m=',i2,' nuse=',i3) 
c END DEBUG
	      if (nuse(m).ge.nuse_limit) then
	         call angstrom(lambda,b_aer,sigb,alfa,sigalfa,m,wrkchan,
     &                         ierr_lang,nb_aer,chi2)
	      endif
	    
	   enddo

*****6******************************************************************
* Print the original and filtered data points, to inspect the 
* quality of the cloud filtering in a plot.
*****6******************************************************************
       call write_filtered_points(xorg,yorg,norg,x,y,nl,
     & fname)

*****6******************************************************************
* Here the results of one day are printed. 
* The main results are the tables of *b_aer* and *alfa* per channel.
* The various results are e.g. the Langley plots.
*****6******************************************************************
	   do m=1,2
	      do i=1,6
	         wrkchan(6)=0
	         if (((wrkchan(i).eq.1).and.(ierr_lang(m,i).lt.99))
     &                 .and.(nuse(m).ge.nuse_limit)) then
	            call write_main(i,index,fname,a,b,
     &                   siga,sigb,sd,b_aer,lambda,alfa,sigalfa,
     &                   proc,m,daynr,yearfrac,ierr_lang,nb_aer,chi2)
	         endif
	      enddo
	   enddo

 1000 enddo
      end

      subroutine auxdata(fname,defozcol,defpsurf,ozcol,psurf)
*****6******************************************************************
* Reads auxiliary data file *readaux.out*.
*
* Produces: ozone column *ozcol* and surface pressure *psurf*.
*
*****6******************************************************************
      implicit none

      integer i,j,j1,j2,dimmax
      parameter(dimmax=5000)
      real ozcol,psurf,defozcol,defpsurf
      real xozone1, xozone2
      real xpres1,xpres2,xrelhum,xwinddir,xwindspeed,xtemp,xvis
      character*10 fname
      character*6 datechar1,datechar2
      character*2 path1
      character*2 dummy 
      logical ex
*****6******************************************************************
      path1='./'

*****6******************************************************************
* Open the aux file.
*****6******************************************************************
      inquire(file=path1//'readaux.out',exist=ex)
      if (.not.ex) then
         write(*,'(3a)') 'File readaux.out does not exist!'
         ozcol=defozcol
         psurf=defpsurf
         return
      else
         open(2,file=path1//'readaux.out')
      endif
*****6******************************************************************
* Read the aux file.
* The ozone and pressure data of AM and PM are averaged,
* because the difference will be small. (??)
*****6******************************************************************
      read(2,9) 
      read(2,9) 
      do i=1,dimmax,2
         read(2,10,end=98) datechar1,j1,xozone1,xpres1,
     &             xtemp,xrelhum,xwindspeed,xwinddir,xvis
         read(2,10,end=98) datechar2,j2,xozone2,xpres2,
     &             xtemp,xrelhum,xwindspeed,xwinddir,xvis
         if ((datechar1.eq.fname(1:6)).and.(datechar2.eq.fname(1:6))) 
     &   then
            if (xozone1.eq.0.) xozone1=defozcol
            if (xozone2.eq.0.) xozone2=defozcol
            if (xpres1.eq.9999.) xpres1=defpsurf
            if (xpres2.eq.9999.) xpres2=defpsurf
            ozcol=0.5*(xozone1+xozone2)
            psurf=0.5*(xpres1+xpres2)
            go to 98
         endif
      enddo
  98  continue
*****6******************************************************************
* Close the aux file.
*****6******************************************************************
      close (2)

      return
*****6******************************************************************
    9 format(a2)
   10 format(a6,1x,i1,7(1x,f6.1))
      end

	subroutine readf(x,y,nl,xorg,yorg,norg,date,dayfr,nmax,
     &  maxsamples,ntot,alln,allx,SPC,SPD,SPA,SPB,fname)
*****6******************************************************************
* Reads the data from the dayfile.
*
*	day in the year,
*	fraction of the day,
*	cosine of the zenith angle,
*	output of the six channels
*
*Produces
*	airmass arrays for morning and afternoon
*	the calibrated measurement of irrandiance
* 
*****6******************************************************************
	implicit none

	integer m,n,i,j,nl(2),norg(2),maxsamples,dimmax,date
	integer nmax,ntot(2),nchan,alln(2)
	parameter(nchan=6)
	parameter(dimmax=5000)
	real x(2,dimmax,nchan),y(2,dimmax,nchan),pi
	real xorg(2,dimmax,nchan),yorg(2,dimmax,nchan)
	real yy(dimmax,nchan),allx(2,dimmax),angle
	double precision newairm
	real costheta(0:dimmax),theta(dimmax),dayfr(2,dimmax)
	real SPC,SPD,SPA(7),SPB(7),sprtime(dimmax),time(dimmax)
	character*10 fname

*****6******************************************************************
* Below the calibration constants, as given by the manufacturer are set.
*****6******************************************************************

      SPA(1) = 2.412917          
      SPB(1) = 0.001200
      SPC =  1000.0
      SPD =     0.0           
        
      SPA(2)= 1.0417198
      SPB(2)= -0.0009700

      SPA(3)=2.60924
      SPB(3)=-0.00040
	
      SPA(4)=0.92835
      SPB(4)=0.00050

      SPA(5)= 1.32079
      SPB(5)= -0.00197

      SPA(6)=9.28147
      SPB(6)=0.00030	
	
      pi=acos(-1.)

   
	do j=1,17
	   read(1,'(a)')
	enddo
	nl(1)=0
	costheta(0)=0.
	do n=1,maxsamples
	   read(1,*,end=99) date,sprtime(n),costheta(n),(yy(n,i),i=1,6)
c	   if (costheta(n).gt.costheta(n-1)) nl(1)=nl(1)+1
	   if (sprtime(n).lt.0.5) nl(1)=nl(1)+1
	   time(n)=sprtime(n)*24*3600
	enddo
  99    continue


*****6******************************************************************
* The data are divided in an AM and a PM part. The number of lines,
* corresponding to AM is nl(1) and the number of lines corresponding
* to PM is nl(2).
* The irradiances are put in arrays y(m,n,i) where m indicates AM or PM,
* *n* counts the data lines and *i* counts the channels.
* The airmasses corresponding to each datapoint are stored
* in arrays x(m,n,i). The airmasses are calculated according to the 
* formula of Kasten and Young, 1989, applied optics, 28,4735-4738
*****6******************************************************************

c	 open(unit=57,file='thetas')
c	 write(57,'(i10,a)') n,fname
c	  write(57,'(f10.8)') pi
	 

	nmax=n-1
	nl(2)=nmax-nl(1)
	ntot(1)=0
	do n=1,nl(1)
	   dayfr(1,n)=sprtime(n)
	   theta(n)=acos(costheta(n))

	 	   
******************************************************************
*
* Convert radians to degrees so that airmass function is valid!
*
*****************************************************************
	   theta(n)= (theta(n)*180.)/pi
	   
*****************************************************************
* The next few lines are only used once in studying the 
* alignment error. An extra zennith angle of half a degree 
* is once added and once substracted.
******************************************************************
c	    theta(n)=theta(n) + 0.5
c	    theta(n)=(theta(n)*pi)/180.
c	    costheta(n)=cos(theta(n))

	   do i=1,6
c              x(1,n,i)=1./(costheta(n)+0.50572*(90.-theta(n)+6.07995)
c     &         **(-1.6364))


	      x(1,n,i)=(1.002432*(costheta(n)**2) + 0.148386*
     &	      costheta(n) + 0.0096467)/(costheta(n)**3 + 
     &	       0.149864*(costheta(n)**2) + 0.0102963*costheta(n)
     &	        + 0.000303978)  	      



c     counter=1.002432*(T**2) + 0.148386*T + 0.0096467
c      numerator= T**3 + 0.149864*(T**2) + 0.0102963*T + 0.000303978
c      newairm=counter/numerator




c	       angle=costheta(n)
c	       x(1,n,i)=newairm(angle)

	      y(1,n,i)=yy(n,i)
	   enddo
	   if(x(1,n,1).gt.2..and.x(1,n,1).lt.6.) then
	      ntot(1)=ntot(1)+1
	      allx(1,ntot(1))=x(1,n,1)
	   endif
	enddo
	


	ntot(2)=0
	do n=1,nl(2)
	   dayfr(2,n)=sprtime(n+nl(1))
	   theta(n+nl(1))=acos(costheta(n+nl(1)))
	   
c	   if (theta(n+nl(1)).gt.((0.5)*pi) then



c	   write(57,'(f12.4)') theta(n)
c	   write(57,907)
c 907       format ('onzin3')

******************************************************************
*
* Convert radians to degrees so that airmass function is valid!
*
*****************************************************************
	   theta(n+nl(1))= (theta(n+nl(1))*180.)/pi
	   
*****************************************************************
* The next few lines are only used once in studying the 
* alignment error. An extra zennith angle of half a degree 
* is once added and once substracted.
******************************************************************
c	    theta(n+nl(1))=theta(n+nl(1)) + 0.5
c	    theta(n+nl(1))=(theta(n+nl(1))*pi)/180.
c	    costheta(n+nl(1))=cos(theta(n+nl(1)))
	   
	   do i=1,6
c	      x(2,n,i)=1.0/(costheta(n+nl(1)) + 0.50572*
c    &        (90.-theta(n+nl(1))+6.07995)**(-1.6364))
     
     
	      x(2,n,i)=(1.002432*(costheta(n+nl(1))**2) + 0.148386*
     &	      costheta(n+nl(1)) + 0.0096467)/(costheta(n+nl(1))
     &	       **3 + 0.149864*(costheta(n+nl(1))**2) + 
     &	       0.0102963*costheta(n+nl(1))+ 0.000303978)  	      

	       
c	       angle=costheta(n)
c	       x(2,n,i)=newairm(angle)
	       
	      y(2,n,i)=yy(n+nl(1),i)
	   enddo
c
c The PM airmass retriction to the maximum value 3.1 
c is replaced by maximum value 6, on 10 Feb. 2003:
c
	   if(x(2,n,1).gt.2..and.x(2,n,1).lt.6.) then
	      ntot(2)=ntot(2)+1
	      allx(2,ntot(2))=x(2,n,1)
	   endif
	enddo

c
c Put the original x and y data in arrays xorg and yorg:
c
	do m=1,2
           norg(m)=nl(m)
	   do n=1,nl(m)
	      do i=1,6
	         xorg(m,n,i)=x(m,n,i)
	         yorg(m,n,i)=y(m,n,i)
	      enddo
	   enddo
	enddo

c
c Calibrate the signals using the factory constants.
c
	do m=1,2
	   do n=1,nl(m)
	      do i=1,6
	         y(m,n,i)=((y(m,n,i))-(SPB(i)*SPC)-SPD)/
     &                         (SPA(i)*SPC)
	      enddo
	   enddo
	enddo

c	write(*,900) fname, date, ntot(1), ntot(2)

c       do m=1,2
c          do n=1,ntot(m)
c             write(22,*)allx(m,n)
c          enddo
c       enddo

  900   format('dayfile: ',a10,' daynumber: ',i3,' data points ',
     &         '(AM, PM):',i4,i4)
	return
	end

	subroutine wrdayvariation(direct,airm,dayfr,nmonth,
     &                            nday,nn,nl,od,fluxfac,fname,m,i)

*****6******************************************************************
* In this subroutine the daily variation of the optical thickness is
* calculated with the calibration constants of the manufacturer.
*****6******************************************************************
*****6******************************************************************
* IMPORTANT:
* BECAUSE OF A 10% ERROR IN THE CALIBRATION, 
* THE ERRORS IN THE CALCULATED
* VALUES OF THE OPTICAL THICKNESS IN THIS WAY ARE VERY LARGE.
* (ABOUT 40 %).
*****6******************************************************************

	implicit none

	integer nchan
	parameter(nchan=6)
	integer nmax,i,n,m,dimmax,nn(2),nl(2,nchan)
	parameter(dimmax=5000)
	integer crit(2,dimmax,nchan)
	real airm(2,dimmax,nchan),dayfr(2,dimmax)
	real xlon,xlat,nmonth,nday
	parameter(xlon=-5.18)
	parameter(xlat=52.10)
	real dist,xmu0,xphi,fluxfac
	real I_ch(6)
	real od(2,dimmax,nchan)
	real direct(2,dimmax,7)
	character*10 fname
	character*18 varfname

	do m=1,2
	   do i=1,6
	      do n=1,nn(m)
	         crit(m,n,i)=1
	         if(direct(m,n,3).lt.0.15.or.direct(m,n,i).le.0.) then
	            crit(m,n,i)=0
	         endif
	      enddo
	   enddo
	enddo

	do m=1,2
	   do i=1,6
	      nl(m,i)=0
	      do n=1,nn(m)
	         if(crit(m,n,i).eq.1)then
		    nl(m,i)=nl(m,i)+1
		    airm(m,nl(m,i),i)=airm(m,n,i)
	            dayfr(m,nl(m,i))=dayfr(m,n)
		    direct(m,nl(m,i),i)=direct(m,n,i)
	         endif
	      enddo
	   enddo
	enddo


       do m=1,2
          do i=1,6
             call dayvariation(direct,airm,nmonth,
     &                    nday,nl,crit,od,fluxfac,m,i)
          enddo
       enddo


c *********************************************************************
c   WRITING THE DAILY VARIATION RESULTS IS COMMENTED TO SAVE SPACE:
c   (P.S., 24 Aug 1998)
c *********************************************************************
c       open(31,file='RESULTS/'//fname(1:6)//'am.var')
c       open(32,file='RESULTS/'//fname(1:6)//'pm.var')
c       do m=1,2
c          do n=1,nl(m,3)
c             write(30+m,22)dayfr(m,n),airm(m,n,3),
c     &        od(m,n,1),od(m,n,3),od(m,n,4),od(m,n,6)
c          enddo
c       enddo
c
c  22   format(f10.5,f10.5,f10.5,f10.5,f10.5,f10.5,f10.5)
*********************************************************************
	return
	end

	subroutine clcdate(daynr,yr,month,day,year,yearfrac)
	implicit none

*****6******************************************************************
* This subroutine calculates the month and day from the daynumber, 
* and the decimal year fraction, *yearfrac*.
*****6******************************************************************

        integer nmonths
        parameter (nmonths=12)
	integer i,j,m,n,yr,month,day,daynr
	integer days(nmonths+1), totday
	real tst,year,yearfrac,fullyear

	data days /0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

        fullyear=365.

*****6******************************************************************
* Test for leap year
* NOTE: 2000 was a leap year !!
*****6******************************************************************
	tst=float(yr)/4
	if ((float(tst)-int(tst)).lt.0.01) then
	   days(3)=29
           fullyear=366.
	endif
	     
*****6******************************************************************
* Calculate. 
*****6******************************************************************
	totday=0
	do n=1,nmonths
	   totday=totday+days(n)
	   if(daynr-totday.le.days(n+1)) goto 20
	enddo
   20   continue
	month=n
	day=daynr-totday

        yearfrac=year + (float(daynr)-0.5)/fullyear

	return
	end

	subroutine firstcheck(x,y,nl,m,i,fname)
	implicit none

*****6******************************************************************
* Here the data points are rejected that do not lie in the airmass
* range 2-6 (see Harrison and Michalsky, 1994).
c
c 10 Feb. 2003: The PM airmass restriction to the maximum value 3.1 
c (which was only relevant to the 1997/1998 data) 
c is now replaced by the maximum value 6.
*****6******************************************************************

	integer dimmax,nchan
	parameter (nchan=6)
	parameter(dimmax=5000)
	integer m,n,i,nl(2,nchan),nn(2,nchan),wrkchan(nchan)
	integer crit(2,dimmax,nchan)
	real x(2,dimmax,nchan), y(2,dimmax,nchan)
	character*30 lanfname
	character*10 fname

	do n=1,nl(m,i)
	   crit(m,n,i)=1
	enddo

	do n=1,nl(1,i)
	  if(x(1,n,i).lt.2..or.x(1,n,i).gt.6.) then
	     crit(1,n,i)=0
	  endif
	enddo

	do n=1,nl(2,i)
	  if(x(2,n,i).lt.2..or.x(2,n,i).gt.6.) then
	     crit(2,n,i)=0
	  endif
	enddo

	nn(m,i)=0
	do n=1,nl(m,i)
	   if(crit(m,n,i).eq.1)then
	      nn(m,i)=nn(m,i)+1
	      x(m,nn(m,i),i)=x(m,n,i)
	      y(m,nn(m,i),i)=y(m,n,i)
	   endif
	enddo
	nl(m,i)=nn(m,i)


	return
	end

	subroutine firstfilter(x,y,nmax,nl,nsum,m,i)
	implicit none

	integer nchan
	parameter(nchan=6)
	integer i,j,k,m,n,nmax,nl(2,nchan),nn(2,nchan)
	integer nsum,dimmax
        real lim_slope
	parameter(dimmax=5000)
	parameter(lim_slope=0.0)
	real x(2,dimmax,nchan),y(2,dimmax,nchan)
	real xav1,xav2,yav1,yav2
	integer crit(2,dimmax,nchan)
	real slope
*****6******************************************************************
* In this subroutine all datapoints are rejected that lie beneath
* the previous data point in the Langley plot(for decreasing airmass). 
* The points cannot correspond to a stable atmosphere because
* Langley plot has a negative slope.
c
c 14 Feb. 2003: the parameter *lim_slope* is introduced to have more
c flexibility in defining the limiting value for the slope.
*****6******************************************************************

	do n=1,nl(m,i)
	   crit(m,n,i)=1	      
	enddo

*****6******************************************************************
c Change on 14 Feb. 2003 (P. Stammes): 
c
c nsum is set to 1, whereas it was in fact equal to 0 in the version 
c used until Feb. 2003, due to a code error.
c 
c The value nsum=1 has been tested to be better than nsum=0, because it
c is more agressive in removing dips from clouds, and better than nsum=2
c (which removed more good points, but even introduced new dips !).
*****6******************************************************************
	nsum=1

	do n=1,nl(m,i)-nsum
	   xav1=0.
	   yav1=0.
	   xav2=0.
	   yav2=0.
	   if(crit(m,n,i).eq.0) goto 10
	   do j=1,nsum
              xav1=xav1+x(m,n+j-1,i)
              yav1=yav1+y(m,n+j-1,i)
              xav2=xav2+x(m,n+j,i)
              yav2=yav2+y(m,n+j,i)
	   enddo
 	   slope=(yav1-yav2)/(xav1-xav2)
	   if(slope.gt.lim_slope) then
	      do j=1,nsum
		 crit(m,n+j-1,i)=0
	      enddo
	      do k=1,nl(m,i)
		 if ((n+k).ge.nmax) goto 10
		 slope=(y(m,n,i)-y(m,n+k,i))/(x(m,n,i)-x(m,n+k,i))
		 if(slope.gt.lim_slope) crit(m,n+k,i)=0
		 if(slope.lt.lim_slope) goto 10
	      enddo
	   endif
   10	enddo

	nn(m,i)=0
	do n=1,nl(m,i)
	   if(crit(m,n,i).eq.1)then
	      nn(m,i)=nn(m,i)+1
	      x(m,nn(m,i),i)=x(m,n,i)
	      y(m,nn(m,i),i)=y(m,n,i)
	   endif
	enddo
	nl(m,i)=nn(m,i)

	return
	end

	subroutine secondfilter(x,y,nl,nsum,nmonth,nday,m,i)
	implicit none

*****6******************************************************************
* In this subroutine an estimate of the slope is made, using the points
* with constant optical thickness. All points are rejected that
* with a slope that differs more than 50 % from the estimated slope.
*****6******************************************************************

	integer dimmax,nmonth,nday,nchan
	parameter(nchan=6)
	parameter(dimmax=5000)
	integer i,j,m,n,nn(2,nchan),nl(2,nchan)
	integer crit(2,dimmax,nchan),nsum
	real x(2,dimmax,nchan),y(2,dimmax,nchan)
	real slope(2,dimmax,nchan)
	real xav1,xav2,yav1,yav2,od(2,dimmax,nchan)
	real odslope(2,dimmax,nchan),sav(2,nchan)
	real nav(2),rsav(2,nchan),fluxfac
	real xx(2,dimmax,nchan),yy(2,dimmax,nchan)

	do n=1,nl(m,i)
	   crit(m,n,i)=1
	enddo
	
	call dayvariation(y,x,nmonth,nday,nl,crit,od,fluxfac,m,i)

	do n=1,nl(m,i)-1
  	   odslope(m,n,i)=abs((od(m,n,i)-od(m,n+1,i))/
     &     (x(m,n,i)-x(m,n+1,i)))
	enddo

	do n=1,nl(m,i)
	   xav1=0.
	   yav1=0.
	   xav2=0.
 	   yav2=0.
	   do j=1,nsum
	      xav1=xav1+x(m,n+j-1,i)
              yav1=yav1+y(m,n+j-1,i)
              xav2=xav2+x(m,n+j,i)
              yav2=yav2+y(m,n+j,i)
	   enddo
	   slope(m,n,i)=(alog(yav1)-alog(yav2))/(xav1-xav2)
	enddo

	nav(m)=0
	sav(m,i)=0
	do n=1,nl(m,i)-1
	   if(abs(odslope(m,n,i)).lt.0.05)then
	      nav(m)=nav(m)+1
	      xx(m,nav(m),i)=x(m,n,i)
	      yy(m,nav(m),i)=y(m,n,i)
	   endif
	enddo

	if (nav(m).eq.0) nl(m,i)=0
	do n=1,nav(m)-1
	   rsav(m,i)=(alog(yy(m,n,i))-alog(yy(m,n+1,i)))/
     &     (xx(m,n,i)-xx(m,n+1,i))
	   sav(m,i)=sav(m,i)+rsav(m,i)
	enddo
	sav(m,i)=sav(m,i)/nav(m)

	do n=1,nl(m,i)
	   if(slope(m,n,i).lt.1.5*sav(m,i).or.
     &	      slope(m,n,i).gt.0.5*sav(m,i)) then
	      crit(m,n,i)=0
	   endif
 	enddo

	nn(m,i)=0
	do n=1,nl(m,i)
	   if(crit(m,n,i).eq.1)then
	      nn(m,i)=nn(m,i)+1
	      x(m,nn(m,i),i)=x(m,n,i)
	      y(m,nn(m,i),i)=y(m,n,i)
	   endif
	enddo
	nl(m,i)=nn(m,i)
   	
   20   format(f10.5,f10.5,f10.5,f10.5,f10.5)

	return
	end

	subroutine thirdfilter(x,y,nl,m,i)
	implicit none

*****6******************************************************************
* In this subroutine all points are rejected that lay more than 1.5
* standard deviation from the regression line.
*****6******************************************************************

	integer nchan
	parameter(nchan=6)
	integer n,m,nn(2,nchan),nl(2,nchan),i,dimmax
	parameter (dimmax=5000)
	integer crit(2,dimmax,nchan)
	real x(2,dimmax,nchan),y(2,dimmax,nchan),a(2),b(2)
	real xx(dimmax),yy(dimmax)
	real siga(2),sigb(2),chi2(2),q(2),sd(2,nchan),ss

	do n=1,nl(m,i)
	   crit(m,n,i)=1
	enddo

	do n=1,nl(m,i)
	    xx(n)=x(m,n,i)
	    yy(n)=alog(y(m,n,i))
	enddo

	call fit(xx,yy,nl(m,i),ss,0.,a(m), b(m),siga(m),sigb(m),
     &           chi2(m), q(m))
     
c	 write(*,*)

	sd(m,i)=sqrt(chi2(m)/float(nl(m,i)-1))

	do n=1,nl(m,i)
	   if(abs(alog(y(m,n,i))-a(m)-x(m,n,i)*b(m)).
     &        gt.1.5*sd(m,i))then
	      crit(m,n,i)=0
	   endif
	enddo

	 nn(m,i)=0
	 do n=1,nl(m,i)
	    if(crit(m,n,i).eq.1)then
	       nn(m,i)=nn(m,i)+1
	       x(m,nn(m,i),i)=x(m,n,i)
	       y(m,nn(m,i),i)=y(m,n,i)
	    endif
	enddo
	nl(m,i)=nn(m,i)

	do n=1,nl(m,i)
	   crit(m,n,i)=1
	enddo
	do n=nl(m,i)+1,dimmax
	   x(m,n,i)=0
	enddo

	return
	end

        subroutine write_filtered_points(xorg,yorg,norg,x,y,nl,fname)
	implicit none

	integer nchan
	parameter(nchan=6)
	integer i,j,k,m,n,dimmax,nmax,nl(2,nchan),norg(2)
	parameter(dimmax=5000)
	real xorg(2,dimmax,nchan),yorg(2,dimmax,nchan)
	real x(2,dimmax,nchan),y(2,dimmax,nchan)
	character*10 fname
	character*16 amname, pmname

c****6******************************************************************
c In this subroutine the original and filtered data points are plotted
c to check in a plot the selection criteria used in the cloud filter
c subroutines: firstcheck, firstfilter, secondfilter, thirdfilter.
c****6******************************************************************
        write (amname,'(a6,a10)') 'AMdata',fname
        write (pmname,'(a6,a10)') 'PMdata',fname

        open (unit=20, file=amname)
        open (unit=21, file=pmname)

	do m=1,2
           write (19+m,1) fname, m
	   do n=1,norg(m)
	      write (19+m,2) (xorg(m,n,j), yorg(m,n,j), 
     &                      x(m,n,j), y(m,n,j), j=1,4)
	   enddo
	enddo

        close (unit=20)
        close (unit=21)

    1  format ('# date:',a10,'    halfday: m=',i1,' (AM=1, PM=2)',
     &  /,'# ',/,'# airm=air mass',
     & 'Irr=irradiance,org=original,fil=filtered,_n=channel number',/, 
     & '# ',/, 
     & '# airm_org_1   Irr_org_1  airm_fil_1  Irr_fil_1  ', 
     & ' airm_org_2   Irr_org_2  airm_fil_2  Irr_fil_2  ',
     & ' airm_org_3   Irr_org_3  airm_fil_3  Irr_fil_3  ',
     & ' airm_org_4   Irr_org_4  airm_fil_4  Irr_fil_4')
c      & '  ao_3 Io_3 af_3 If_3 ao_4 Io_4 af_4 If_4', 
c      & '  ao_5 Io_5 af_5 If_5 ao_6 Io_6 af_6 If_6') 
    2  format (4(1x,e11.4,1x,e11.4,1x,e11.4,1x,e11.4))

        return
        end

	subroutine langley(x,y,nl,fname,daynr,ra,rsiga,rb,rsigb,
     &  b_aer,sd,I_ex,fluxfac,b_ozon,err_oz,bral,m,i,toterr,ierr_lang)
c
c Perform Langley regression to obtain optical thickness and
c extraterrestrial solar irradiance.
c
c Error codes are generated in case the fit yields a positive slope
c (thus, negative b_aer).
c (P.S., 2 Aug. 1999)
c
        implicit none

	integer nchan
	parameter (nchan=6)
	integer i,j,m,n,nl(2,nchan),nn(2,nchan),dimmax,nhead
	integer daynr,pdef,odef
	integer wrkchan(nchan)
	parameter(nhead=12)
	parameter (dimmax=5000)
	integer ierr_lang(2,nchan)
	real x(2,dimmax,nchan),y(2,dimmax,nchan)
	real rsiga(2,nchan),rsigb(2,nchan)
	real xx(dimmax),yy(dimmax),siga(2),sigb(2)
	real sd(2,nchan),proc(2,nchan)
	real chi2(2),q(2),a(2),b(2),ra(2,nchan),rb(2,nchan),ss
	real bral(nchan),err_oz(nchan)
	real alf(2),b_ozon(nchan),b_aer(2,nchan),fluxfac
        real toterr(2,nchan)
	real rb_aer(nchan),I_ex(2,nchan)
	real SPA(nchan),SPB(nchan)
	character*10 fname
	character*120 header(nhead)
	character*30 taufname
	character*30 resfname(nchan) 
	logical ex

	SPA(1) = 2412.917          
	SPB(1) = 1.200
       
	SPA(2)= 1041.7198
	SPB(2)= -0.9700

	SPA(3)=2609.24
	SPB(3)=-0.40
	
	SPA(4)=928.35
	SPB(4)=0.50

	SPA(5)= 1320.79
	SPB(5)= -0.197

	SPA(6)=9281.47
	SPB(6)=0.30

	do n=1,nl(m,i)
	   y(m,n,i)=y(m,n,i)*SPA(i)+SPB(i)
	enddo

	do n=1,nl(m,i)
	   xx(n)=x(m,n,i)
	   yy(n)=alog(y(m,n,i))	
	enddo

	call fit(xx,yy,nl(m,i),ss,0.,a(m),b(m),siga(m),sigb(m),
     &           chi2(m),q(m))

	ra(m,i)=a(m)
	rsiga(m,i)=siga(m)
	rb(m,i)=-b(m)

	rsigb(m,i)=sigb(m)
	toterr(m,i)=sqrt(rsigb(m,i)*rsigb(m,i)+err_oz(i)*err_oz(i))

	sd(m,i)=sqrt(chi2(m)/float(nl(m,i)-1))

	b_aer(m,i)=rb(m,i)-bral(i)-b_ozon(i)
c
c If the retrieved aerosol optical thickness b_aer is negative,
c then make b_aer equal to a small positive value, 0.00001.
c If in this case the possible error on b_aer is larger than |b_aer|, 
c ierr_lang=1, else ierr_lang=2.
c
        if (b_aer(m,i).lt.0.0) then
           if (toterr(m,i).gt.abs(b_aer(m,i))) then
              ierr_lang(m,i)=1
           else
              ierr_lang(m,i)=2
           endif
           b_aer(m,i)=0.00001
        else
           ierr_lang(m,i)=0
        endif

	I_ex(m,i)=(exp(ra(m,i)))/fluxfac
	ra(m,i)=alog(I_ex(m,i))
c	write(*,*) 'endfit',ra(m,i),rb(m,i),b_aer(m,i)

	return
	end

	subroutine dayvariation(direct,airm,nmonth,
     &                          nday,nl,crit,od,fluxfac,m,i)

	implicit none

	integer nchan
	parameter (nchan=6)
	integer nmax,i,n,m,dimmax,nn(2,nchan),nl(2,nchan)
	parameter(dimmax=5000)
	integer crit(2,dimmax,nchan)
	real channel(dimmax,nchan),airm(2,dimmax,nchan)
	real xlon,xlat,nmonth,nday
	parameter(xlon=-5.18)
	parameter(xlat=52.10)
	real dist,xmu0,xphi,fluxfac
	real I_ch(6)
	real x(2,dimmax,nchan),y(2,dimmax,nchan)
	real od(2,dimmax,nchan)
	real direct(2,dimmax,7)
	character*10 fname

	I_ch(1)=1.0959747
	I_ch(2)=1.8174066
	I_ch(3)=1.5166416*1.106
	I_ch(4)=1.1897208*1.102
	I_ch(5)=0.9844523
	I_ch(6)=0.8059377

c	I_ch(1)=2645.695985
c	I_ch(2)=1892.25844
c	I_ch(3)=4376.353813
c	I_ch(4)=1941.422898
c	I_ch(5)=1298.284753
c	I_ch(6)=7480.586584

	call sunpos(xlon,xlat,nmonth,nday,0,0,0,
     &  dist,xmu0,xphi,fluxfac)

	I_ch(i)=alog(I_ch(i)*fluxfac)

 	do n=1,nl(m,i)
	   od(m,n,i)=abs((I_ch(i)-(alog(direct(m,n,i))))/
     &     airm(m,n,i))
	enddo

	return
	end

	subroutine angstrom(lambda,b_aer,sigb,alfa,sigalfa,m,wrkchan,
     &                         ierr_lang,nb_aer,chi2)
*****6******************************************************************
c Here the Angstrom exponent, alfa, is determined by fitting b_aer to an
c exponential function of wavelength: 
c
c     b_aer (lambda) = beta x (lambda)^-alfa
c
c This formula is rewritten as follows, to be able to use 
c the linear-regression routine FIT from Press et al.:
c
c <=> log b_aer = log(beta) - alfa log(lambda)
c
c The error (sigma) in b_aer, sigb, is used as follows in the
c linear regression:
c
c     log (b_aer + sigb) = log (b_aer) + log (1 + sigb/b_aer)
c
c So, the logarithmic error is: rsigb_aer = log (1 + sigb/b_aer)
*****6******************************************************************
	implicit none

	integer nchan
	parameter (nchan=6)
	integer i,j,k,m,n,nb_aer(2),wrkchan(nchan),ierr_lang(2,nchan)
	real b_aer(2,nchan),rb_aer(nchan)
	real sigb(2,nchan)
	real lambda(nchan),rlambda(nchan),alfa(2),beta(2)
	real rsigb_aer(nchan)
	real siga(2),chi2(2),q(2),ss,a(2), sigalfa(2)

*****6******************************************************************
* Initialize the temporary arrays.
*****6******************************************************************
        do i=1,nchan
           rlambda(i)=0.0
           rb_aer(i)=0.0
           rsigb_aer(i)=0.0
        enddo
           
*****6******************************************************************
* The watervapour channel is not used for determining *alfa*, so
* wrkchan(6)=0
*****6******************************************************************
	wrkchan(6)=0

*****6******************************************************************
* Prepare the x- and y-arrays for the linear regression.
*****6******************************************************************
	n=0
        do i=1,nchan
           if ((wrkchan(i).eq.1).and.(ierr_lang(m,i).le.1)) then
	      n=n+1
              rb_aer(n)=alog(b_aer(m,i))
              rsigb_aer(n)=alog(1+(sigb(m,i)/b_aer(m,i)))
	      rlambda(n)=-alog(lambda(i)/1.E09)
           endif
        enddo
        nb_aer(m)=n

        if (nb_aer(m).le.2) then
c           write (*,'(a)') 'nb_aer too small !!'
           alfa(m)=99.
           sigalfa(m)=99.
           chi2(m)=-99.
           beta(m)=-99.
        
           return
        endif

*****6******************************************************************
* The standard deviations of b_aer are included in the regression 
* (MWT=1).
*****6******************************************************************
        call fit(rlambda,rb_aer,n,rsigb_aer,1.,a(m),alfa(m),
     &           siga(m),sigalfa(m),chi2(m),q(m))
        beta(m)=exp(a(m))

c DEBUG:
c        write (*,900) m
c        write (*,901) (rlambda(k),k=1,n)
c        write (*,902) (rb_aer(k),k=1,n)
c        write (*,903) (rsigb_aer(k),k=1,n)
c  900 format('m:',i2) 
c  901 format('rlambda:',5(1x,f6.1)) 
c  902 format('rb_aer:',5(1x,f7.4)) 
c  903 format('rsigb_aer:',5(1x,f7.4)) 
c        write (*,904) alfa(m), sigalfa(m), chi2(m) 
c  904 format('alfa: ',f7.4,' sigalfa: ',f7.4,' chi2: ',e9.2) 
c
c END DEBUG:


*****6******************************************************************
* wrkchan(6) must again be set to 1
*****6******************************************************************
	wrkchan(6)=1

	return
	end

	subroutine write_main(ic,index,fname,a,b,
     &     siga,sigb,sd,b_aer,lambda,alfa,sigalfa,
     &     proc,m,daynr,yearfrac,ierr_lang,nb_aer,chi2)

*****6******************************************************************
* This subroutine prints the main results per channel to 
* files "channelX.out", where X is the channel number *ic*.
*****6******************************************************************

	implicit none

	integer nchan,nhead, index
	parameter(nhead=14)
	parameter (nchan=6)
	integer i,j,m,n,ic,nb_aer(2)
	integer ierr_lang(2,nchan)
	integer daynr,refch
	real a(2,nchan),b(2,nchan),siga(2,nchan),sigb(2,nchan)
	real b_aer(2,nchan),alfa(2),chi2(2),sigalfa(2)
	real sd(2,nchan),lambda(6)
	real proc(2,nchan)
        real yearfrac
	character*10 fname
	character*120 header(nhead)
	character*30 filename(nchan)
	logical ex

*****6******************************************************************
* Printing the header of the file.
*****6******************************************************************
	write(header(2),'(a13,i2,a30,a26)')
     &  '# file index=',index,' (1: O3 + p from readaux.out; 0: O3 + p',
     &  ' values read from spuv.in)'
	header(3)='# m=1 / 2: morning / afternoon'
	header(4)='# % : percentage of points useful for regression'
	header(5)='# sd: standard deviation of Langley regression'
	header(6)='# I(0) = extraterr. solar irradiance in counts'
        header(7)='# sigX = standard deviation of X '
	header(8)='# b = optical thickness (tot=total, aer=aerosol)'
	header(9)='# nb_aer = number of wavelengths to derive alpha'
	header(10)='# alpha = Angstrom coefficient'
	header(11)='# chi2  = chi-square of Angstrom relation fit'
	header(12)='# error (in Langley plot): 0=none,1=negative b_aer,
     &   within error bar, 2=negative b_aer, outside error bar'
	header(13)='# '
	write(header(14),'("# date m day yearfrac  %    sd'//
     &  '     lnI(0)  siglnI0   b_tot  sigb_tot '//
     &  ' b_aer nb_aer  alfa   sigalfa   chi2   error")')  

*****6******************************************************************
* Printing the results.
*****6******************************************************************
	write(header(1),'(a12,f6.1,a3)')
     &  '# wavelength=',lambda(ic),' nm'
	write(filename(ic),'(a8,a7,i1,a4)')
     &        'RESULTS/','channel',ic,'.out'
	inquire(file=filename(ic), exist=ex)
	if (.not.ex) then
	   open(unit=60+ic, file=filename(ic),
     &        access='append')
	   do j=1,nhead
	      write(60+ic,'(a)') header(j)
	   enddo
	else
	   open(unit=60+ic, file=filename(ic),
     &              access='append')
	endif
	write(60+ic,20) fname(1:6),m,daynr,yearfrac,proc(m,ic),sd(m,ic),
     &        a(m,ic),siga(m,ic), b(m,ic),sigb(m,ic),
     &        b_aer(m,ic),nb_aer(m),alfa(m),sigalfa(m),chi2(m),
     &        ierr_lang(m,ic)

*****6******************************************************************
   20   format(a6,1x,i1,1x,i3,1x,f8.3,1x,f4.1,1x,f6.4,1x,f7.4,1x,
     &         f7.4,1x,f8.5,1x,f7.5,1x,f8.5,1x,i3,3x,f7.4,1x,
     &         f7.4,1x,e9.2,1x,i3)
	return
	end

	subroutine sunpos(xlon,xlat,nmonth,nday,nhour,nmin,nsec,
     &  dist,xmu0,xphi,fluxfac)
c-----------------------------------------------------------------------
c Given the latitude and longitude of the location on Earth (in degrees)
c and the date and time (in GMT) this module computes
c dist - the relative Sun-Earth distance
c xmu0 - the cosine of the solar zenith angle
c xphi - the solar azimuth angle
c according to A.C. Velds, "Zonnestraling in Nederland", 1992
c-----------------------------------------------------------------------
 
      dimension nm(13)
      data nm /0,31,59,90,120,151,181,212,243,273,304,334,365/

      pi=4.*atan(1.)
      rtd=180./pi
      dtr=1./rtd
      daynum=float(nm(nmonth)+nday)
c
c relative Sun-Earth distance, p.99
c
      eta=2.*pi/365.
      fluxfac=1.000110
     1  +0.034221*cos(eta*daynum)+0.000719*cos(2.*eta*daynum)
     2  +0.001280*sin(eta*daynum)+0.000077*sin(2.*eta*daynum)
      dist=1./sqrt(fluxfac)
c
c solar declination
c
      ed=2.*pi*daynum/365.
      delta=0.006918
     1  -0.399912*cos(ed)-0.006758*cos(2.*ed)-0.002697*cos(3.*ed)
     2  +0.070257*sin(ed)+0.000907*sin(2.*ed)+0.001480*sin(3.*ed) 
c
c equation of time
c
      et=2.*pi*daynum/366.
      eq=0.0072*cos(et)-0.0528*cos(2.*et)-0.0012*cos(3.*et)
     1  -0.1229*sin(et)-0.1565*sin(2.*et)-0.0041*sin(3.*et) 
c
c solar position p.132
c
  
      time=float(nhour)+float(nmin)/60.+float(nsec)/3600.
      omega=(360./24.)*(time-(xlon/15.)+eq-12.)*dtr
      sinh=sin(delta)*sin(xlat*dtr)+cos(delta)*cos(xlat*dtr)*cos(omega)
      solel=asin(sinh)
      xmu0=cos(pi/2.-solel)
      if ((abs(solel-pi/2.).gt.1.e-5).and.
     1  (abs(solel-3.*pi/2.).gt.1.e-5)) then
        cospsi=(-sin(delta)*cos(xlat*dtr)+
     1    cos(delta)*sin(xlat*dtr)*cos(omega))/cos(solel)
        sinpsi=cos(delta)*sin(omega)/cos(solel)
        xphi=rtd*atan2(sinpsi,cospsi)
      endif

      return
      end


 
      SUBROUTINE FIT(X,Y,NDATA,SIG,MWT,A,B,SIGA,SIGB,CHI2,Q)            
c
c See Press et al., 1986, page 508.
c
C Given a set of NDATA points X(I), Y(I), with standard deviations
C SIG(I), fit them to a straight line y=a+bx by minimizing chi-square.
C Returned are A, B, and their respective probable uncertainties
C SIGA and SIGB, the chi-square CHI2 and the goodness of fit probability
C Q (that the fit would have chi-square this large or larger).
C If MWT=0 on input, then the standard deviations are assumed to be
C unavailable: Q is returned as 1.0 and the normalization of CHI2
C is to unit standard deviation on all points.
C
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA)                            
      SX=0.                                                             
      SY=0.                                                             
      ST2=0.                                                            
      B=0.                                                              
      IF(MWT.NE.0) THEN                                                 
        SS=0.                                                           
        DO 11 I=1,NDATA                                                 
          WT=1./(SIG(I)**2)                                             
          SS=SS+WT                                                      
          SX=SX+X(I)*WT                                                 
          SY=SY+Y(I)*WT                                                 
11      CONTINUE                                                        
      ELSE                                                              
        DO 12 I=1,NDATA                                                 
          SX=SX+X(I)                                                    
          SY=SY+Y(I)                                                    
12      CONTINUE                                                        
        SS=FLOAT(NDATA)                                                 
      ENDIF                                                             
      SXOSS=SX/SS                                                       
      IF(MWT.NE.0) THEN                                                 
        DO 13 I=1,NDATA                                                 
          T=(X(I)-SXOSS)/SIG(I)                                         
          ST2=ST2+T*T                                                   
          B=B+T*Y(I)/SIG(I)                                             
13      CONTINUE                                                        
      ELSE                                                              
        DO 14 I=1,NDATA                                                 
          T=X(I)-SXOSS                                                  
          ST2=ST2+T*T                                                   
          B=B+T*Y(I)                                                    
14      CONTINUE                                                        
      ENDIF                                                             
      B=B/ST2                                                           
      A=(SY-SX*B)/SS                                                    
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)                                 
      SIGB=SQRT(1./ST2)                                                 
      CHI2=0.                                                           
      IF(MWT.EQ.0) THEN                                                 
        DO 15 I=1,NDATA                                                 
          CHI2=CHI2+(Y(I)-A-B*X(I))**2                                  
15      CONTINUE                                                        
        Q=1.                                                            
        SIGDAT=SQRT(CHI2/(NDATA-2))                                     
        SIGA=SIGA*SIGDAT                                                
        SIGB=SIGB*SIGDAT                                                
      ELSE                                                              
        DO 16 I=1,NDATA                                                 
          CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2                         
16      CONTINUE                                                        
        Q=GAMMQ(0.5*(NDATA-2),0.5*CHI2)                                 
      ENDIF                                                             
      RETURN                                                            
      END                                                               

      FUNCTION GAMMQ(A,X)                                               
      IF(X.LT.0..OR.A.LE.0.)PAUSE                                       
      IF(X.LT.A+1.)THEN                                                 
        CALL GSER(GAMSER,A,X,GLN)                                       
        GAMMQ=1.-GAMSER                                                 
      ELSE                                                              
        CALL GCF(GAMMQ,A,X,GLN)                                         
      ENDIF                                                             
      RETURN                                                            
      END                                                               

      SUBROUTINE GSER(GAMSER,A,X,GLN)                                   
      PARAMETER (ITMAX=100,EPS=3.E-5)                                   
c
c  On 30/7/99 eps has been changed from 3.e-7 to 3.e-5. (P.S.)
c
      GLN=GAMMLN(A)                                                     
      IF(X.LE.0.)THEN                                                   
        IF(X.LT.0.)PAUSE                                                
        GAMSER=0.                                                       
        RETURN                                                          
      ENDIF                                                             
      AP=A                                                              
      SUM=1./A                                                          
      DEL=SUM                                                           
      DO 11 N=1,ITMAX                                                   
        AP=AP+1.                                                        
        DEL=DEL*X/AP                                                    
        SUM=SUM+DEL                                                     
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1                             
11    CONTINUE                                                          
      PAUSE 'A too large, ITMAX too small'                              
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)                                   
      RETURN                                                            
      END                                                               

      SUBROUTINE GCF(GAMMCF,A,X,GLN)                                    
      PARAMETER (ITMAX=100,EPS=3.E-7)                                   
      GLN=GAMMLN(A)                                                     
      GOLD=0.                                                           
      A0=1.                                                             
      A1=X                                                              
      B0=0.                                                             
      B1=1.                                                             
      FAC=1.                                                            
      DO 11 N=1,ITMAX                                                   
        AN=FLOAT(N)                                                     
        ANA=AN-A                                                        
        A0=(A1+A0*ANA)*FAC                                              
        B0=(B1+B0*ANA)*FAC                                              
        ANF=AN*FAC                                                      
        A1=X*A0+ANF*A1                                                  
        B1=X*B0+ANF*B1                                                  
        IF(A1.NE.0.)THEN                                                
          FAC=1./A1                                                     
          G=B1*FAC                                                      
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1                             
          GOLD=G                                                        
        ENDIF                                                           
11    CONTINUE                                                          
      PAUSE 'A too large, ITMAX too small'                              
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G                                    
      RETURN                                                            
      END          

      FUNCTION GAMMLN(XX)                                               
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER                          
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,          
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/     
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/                              
      X=XX-ONE                                                          
      TMP=X+FPF                                                         
      TMP=(X+HALF)*LOG(TMP)-TMP                                         
      SER=ONE                                                           
      DO 11 J=1,6                                                       
        X=X+ONE                                                         
        SER=SER+COF(J)/X                                                
11    CONTINUE                                                          
      GAMMLN=TMP+LOG(STP*SER)                                           
      RETURN                                                            
      END  
                                                                   
      FUNCTION GAMMP(A,X)                                               
      IF(X.LT.0..OR.A.LE.0.)PAUSE                                       
      IF(X.LT.A+1.)THEN                                                 
        CALL GSER(GAMMP,A,X,GLN)                                        
      ELSE                                                              
        CALL GCF(GAMMCF,A,X,GLN)                                        
        GAMMP=1.-GAMMCF                                                 
      ENDIF                                                             
      RETURN
                                                         
      END     
      
      
      FUNCTION newairm(T)
      implicit none
      double precision counter, numerator,newairm
      double precision T
      
      write(*,*) 'new airmass function is used'
      counter=1.002432*(T**2) + 0.148386*T + 0.0096467
      numerator= T**3 + 0.149864*(T**2) + 0.0102963*T + 0.000303978
      newairm=counter/numerator
      return
      end
