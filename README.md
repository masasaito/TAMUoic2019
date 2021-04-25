-------------- READ-ME file for the TAMUoic2019 database ---------------

 This is a database of the single-scattering properties of oriented 
 ice crystals (OIC). The two major oriented particles, horizontally 
 oriented plates (HOP) and columns (HOC), are considered in the 
 database. The organization of this database is described in the 
 ReadMe file.


 Lead developer: Dr. Masanori Saito (masa.saito@tamu.edu)
 PI: 		 Dr. Ping Yang 


 History
 (MM/DD/YYYY)	(Author)	(Contents)
 02/01/2018	Masa Saito	Project initiated.
 12/10/2018 	Masa Saito	HOP completed. 
 03/01/2019	Masa Saito	HOC completed. 
 07/10/2019	Masa Saito	Version 1.0.0 Released.
 04/18/2021	Masa Saito	Encoder and kernel technique 
                                implemented.
 04/24/2021	Masa Saito	Version 1.1.0 Released.


 
1. General information

   The TAMUoic2019 database includes the single-scattering properties
   of oriented ice crystals. This database is useful for users who 
   would like to consider horizontally oriented ice crystals in ice 
   clouds in various radiative transfer simulations such as lidar 
   signals and brightness distributions on the sky. Several example 
   simulations are performed and summarized in the following paper.

   Reference:

     1). Saito, M., and P. Yang (2019), Oriented ice crystals: A 
         single-scattering property database for applications to lidar
         and optical phenomenon simulations, J. Atmos. Sci., 76, 
         2635--2652.



   1.1 Code and Database

     The TAMUoic2019 database consists of the following directories:
     
       src		:: Source code
       src/hparx  	:: HPARX library used for the TAMUoic2019 
                           database encoding
       bin		:: Executable command path
       examples		:: Example namelist files
       examples/output  :: Directory for the output from the test run 
       params           :: Directory for the data-kernel parameter 
                           files

     In addition to these three directories, users will need to obtain 
     the following directory data from zenodo.org (for more information,
     please visit https://sites.google.com/site/masanorisaitophd/
     data-and-resources/tamuoic2019).

       tables	 	:: Data kernels



   1.2 Variable Definition

      W     :: wavelength (micron, um)
      D     :: maximum dimension (um), defined as an average of
               circumscribed diameters of individual particles.
      th_i  :: Incident polar angle (deg) 
      V     :: volume of a particle (um^3)  
      Aproj :: projected area in case of random orientation (um^2)
      Ageo  :: projected area in case of a specific orientation (um^2)
      Qe    :: extinction efficiency
      SSA   :: single-scattering albedo
      g     :: asymmetry factor
      P??   :: phase matrix elements (?? includes 11, 12, 22, 33, 43, 
               and 44)

      Whole steradian integral of P11 is 4*pi:
    
           pi     2*pi
        int    int     P11 d(zenith angle)d(azimuth angle)=4*pi
           0      0

      Other independent phase matrix elements P12, P13, ..., and P44 
      are normalized by P11.



 2. Installation

   First, to install the database creation codes, edit the following 
   files:

     ./src/Common.mk
     ./src/hparx/Common.mk  
     
   In these files, users can specify their compiler (FC; e.g., gfortran)
   and associate options (FCFLAGS).

   Then, command as follows:

     $ cd ./src/hparx
     $ make
     $ cd ../
     $ make install
 
   In the end, the following executable files will be created in the 
   "bin" directory:

     tamuoic2019encode	:: Encode the single-scattering property 
  			   database of oriented ice crystals with the
                           ASCII format.


      
 3. How to run the codes

   3.1 Examples

     Go to the example directory (examples), and users will find 
     several files and directories as follows:
 
       output/ 			  :: a directory for output 
       TAMUoic2019encode_exp1.nml :: an example namelist (single case)
       TAMUoic2019encode_exp2.nml :: an example namelist (size-resolved
                                     database)
       TAMUoic2019encode_exp3.nml :: an example namelist (size and 
                                     incident angle resolved database)

     Simply, type the following command

       $ ../bin/tamuoic2019encode TAMUoic2019encode_exp?.nml
 
     Then, the ASCII format of the user-specific database is stored in
     the "output" directory (by default).  



   3.2 Namelist

     The namelist contains 

     ******* Begin *******
     &TAMUoic2019

        TAMUoic2019path = '../tables' ! (FIXED)
        outpath         = './output'  ! (USER-SPECIFIED)

        nsiz = 1                      ! (USER-SPECIFIED)
        nwav = 1                      ! (USER-SPECIFIED)
        niang= 1                      ! (USER-SPECIFIED)
        siz(1:1) = 1.0000000E+02      ! (USER-SPECIFIED)
        wav(1:1) = 3.5500000E-01      ! (USER-SPECIFIED)
        iang(1:1)= 3.0000000E-01      ! (USER-SPECIFIED)
        hab    = 'OrientedPlate'      ! (USER-SPECIFIED) 
        pol    = 'Vector'             ! (USER-SPECIFIED) 
     /
     *******  End  *******

     Users do not have to change the table configurations marked as 
     (FIXED).

     For other parameters marked as (USER-SPECIFIED) is what users can 
     specify  as explained in the subsections below, including 
     wavelengths, incident polar angles, and maximum diameters. 

     To specify multiple bins in user-specified parameters, a user 
     needs to modify the max. # of elements of a parameter and 
     associate bin parameter. For example, a namelist that produce the 
     database with multiple maximum diameters should have the following
     lines:
 
        nsiz = {USER-SPECIFIED_NUMBER_N}
        siz(1:{USER-SPECIFIED_NUMBER_N}) = {Bin1}, {Bin2}, ..., {Bin_N}

     'hab' indicates the ice crystal habit, which should be either
     "OrientedPlate" for HOP or "OrientedColumn" for HOC.
     'pol' indicates the polarization in the phase matrix, which should 
     be either "Scaler" or "Vector". 

     To run tamuoic2019encode, simply type:

       $ ../bin/tamuoic2019encode {NAMELIST}


 
4. Output 

   Output is 17 ASCII files including:

     isca.dat	:: Geometric and scattering properties 
     P11.dat  :: Scattering phase matrix element P11
     P12.dat  :: Normalized scattering phase matrix element P12/P11
     P13.dat  :: Normalized scattering phase matrix element P13/P11
     P14.dat  :: Normalized scattering phase matrix element P14/P11
     P21.dat  :: Normalized scattering phase matrix element P21/P11
     P22.dat  :: Normalized scattering phase matrix element P22/P11
     P23.dat  :: Normalized scattering phase matrix element P23/P11
     P24.dat  :: Normalized scattering phase matrix element P24/P11
     P31.dat  :: Normalized scattering phase matrix element P31/P11
     P32.dat  :: Normalized scattering phase matrix element P32/P11
     P33.dat  :: Normalized scattering phase matrix element P33/P11
     P34.dat  :: Normalized scattering phase matrix element P34/P11
     P41.dat  :: Normalized scattering phase matrix element P41/P11
     P42.dat  :: Normalized scattering phase matrix element P42/P11
     P43.dat  :: Normalized scattering phase matrix element P43/P11
     P44.dat  :: Normalized scattering phase matrix element P44/P11.


   
   4.1 "isca.dat"

     'isca.dat' consists of nwav*niang*nsiz lines and 9 columns 
     as follows:

     1). wavelength (µm) 
     2). maximum dimension of the particle (µm)
     3). incident polar angle (deg)
     4). volume of particle (um^3) 
     5). geometric projected area (um^2)
     6). projected area for random orientation (um^2)
     7). extinction efficiency
     8). single-scattering albedo 
     9). asymmetry factor

     For example, 'isca.dat' with nwav=2, niang=10, and nsiz=20 consists
     of 2*10*20=400 lines in which the first 20 lines are for various 
     sizes with a first incident polar angle and wavelength. The next 20 
     lines are for various sizes with a second incident polar angle and 
     the first wavelength. The loop structure is as follows:

       do iwav = 1, nwav        ! wavelength
          do iiang = 1, niang   ! incident polar angle 
             do isiz = 1, nsiz  ! maximum diameter
                (9 columns; W, D, th_i, V, Ageo, Aproj, Qext, SSA, g)
             end do
          end do
       end do



   4.2 "P??.dat"

     'P??.dat' consists of nwav*niang*nsiz*naang+1 lines and 498 columns
     corresponding to the scattering polar angle. The first line shows 
     scattering angles (0-180 degree) and the rest of lines shows the  
     corresponding two-dimensional phase matrix function (P11) or their
     normalized counterparts (e.g., P??/P11). The loop structure is as
     follows.

       do iwav = 1, nwav        ! wavelength
          do iiang = 1, niang   ! incident polar angle 
             do isiz = 1, nsiz  ! maximum diameter
                *** 2D phase matrix element block BEGIN ***
                do iaang = 1, naang      ! scattering azimuthal angle            
                   (498 columns; P?? for each scatter polar angle)
                end do
                *** 2D phase matrix element block END ***
             end do
          end do
       end do

  

  5. Requirements

    The lead developer, Dr. Masanori Saito encourages the research
    community to utilize the TAMUoic2019 database for ice cloud 
    studies. The only requirement in regards to utilizing the 
    scattering property database for research is to acknowledge our 
    contribution in papers/presentations by:

      1. cite Saito and Yang (2019) in a relevant section of the 
         main text.
      2. add the following description in Acknowledgement section or 
         Data Availability section:
      
    "The scattering properties are obtained from TAMUoic2019."



  6. Contact Information

     Dr. Masanori Saito (masa.saito@tamu.edu) 
     Department of Atmospheric Sciences, Texas A&M University. 

 
