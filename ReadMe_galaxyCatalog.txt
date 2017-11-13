I/221               The NED Nearby Galaxy Catalog - NNGC (French+ 2018)
================================================================================
The NED Nearby Galaxy Catalog - NNGC - the rest is an example:
     Tucholke H.-J., de Boer K.S., Seitter W.C.
    <Astron. Astrophys. Suppl. Ser., 119, 91-98 (1996)>
    <The Messenger 81, 20 (1995)>
    =1996A&AS..119...91T 
    =1995Msngr..81...20D
================================================================================
ADC_Keywords: Galaxies, nearby ; Combined data

Description:
    The Magellanic Catalogue of Stars (MACS) is based on scans of ESO
    Schmidt plates and contains about 244,000 stars covering large areas
    around the LMC and the SMC. The limiting magnitude is B<16.5m and the
    positional accuracy is better than 0.5" for 99% of the stars. The
    stars of this catalogue were screened interactively to ascertain that
    they are undisturbed by close neighbours.


File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl    Records    Explanations
--------------------------------------------------------------------------------
ReadMe          80          .    This file
lmc.dat         52     175779    The Large Magellanic Cloud
smc.dat         52      67782    The Small Magellanic Cloud
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lmc.dat smc.dat
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1- 27  A27	 ---	 Name	  Preferred name
  29- 56  A27	 ---	 NEDName  NED preferred name
  58- 66  F8.6   ---     z        Redshift
  68- 

   1- 12  A12    ---     MACS     Designation
  14- 15  I2     h       RAh      Right Ascension J2000 , Epoch 1989.0 (hours)
  17- 18  I2     min     RAm      Right Ascension J2000 (minutes)
  20- 25  F6.3   s       RAs      Right Ascension J2000 (seconds)
      27  A1     ---     DE-      Declination J2000 (sign)
  28- 29  I2     deg     DEd      Declination J2000 , Epoch 1989.0 (degrees)
  31- 32  I2     arcmin  DEm      Declination J2000 (minutes)
  34- 38  F5.2   arcsec  DEs      Declination J2000 (seconds)
      40  I1     ---     Npos     Number of positions used
  42- 46  F5.2   mag     Mag      ?=99.00 Instrumental Magnitude
                                        (to be used only in a relative sense)
      48  I1     ---     PosFlag  [0,1] Position Flag   (0: ok,
                                        1: internal error larger than 0.5")
      50  I1     ---     MagFlag  [0,1] Megnitude Flag  (0: ok,
                                        1: bad photometry or possible variable)
      52  I1     ---  BochumFlag *[0] Bochum Flag
--------------------------------------------------------------------------------
Note on BochumFlag: 1 if in Bochum catalog of astrophysical information
    on bright LMC stars) (yet empty)
--------------------------------------------------------------------------------

Author's address:
    Hans-Joachim Tucholke    <tucholke@astro.uni-bonn.de>

================================================================================
(End)            Hans-Joachim Tucholke [Univ. Bonn]                  20-Nov-1995