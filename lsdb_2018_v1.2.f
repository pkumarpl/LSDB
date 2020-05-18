cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 
c  PROGRAM  L S D B  -  AUTOMATIC SETUP OF THE LOCAL COORDINATE SYSTEM
c 
c Program reads PLATON, SHELX INS file, XD MASTER and PARAMETER files, or
c AIM2TAB output file and tries to figure out the local coordinate system
c for each atom based on connectivity and geometrical parameters. The 
c defined local axes can be visualized in PLATON (at least the UNIX 
c version of PLATON).
c 
c In addition it provides an interface to COPPENS LAB pseudoatom databank
c so that pseudoatom parameters can automatically be assigned to atoms in 
c the structure.
c
c Current version has also the ability to automatically extend positions
c of hydrogen atoms along the predetermined vectors to standard/neutron 
c distances according to data published in International Tables
c 
c It can also be used as advanced converter from SHELX to XD as it now
c correctly takes care of FVAR/AFIX/HFIX etc instructions
c
c Unauthorized use of this program is strictly prohibited and subject for
c prosecution under the penalty of laws of the Russian Federation and USA
c 
c Written by Dr. Anatoliy Volkov w/contrib. from M. Messerschmidt,
c K.N. Jarzembska, P. Kumar, P. Dominiak (c) University at Buffalo, 2001-2006,
c Warsaw University, 2006-2018
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Modifications 
c -------------
c 121610: Added exceptions in DetExcept for:
c         Fluorine: PH-F
c         Chlorine: PH-Cl
c         New atom keys added: number of planar rings to which atom belongs,
c         RING, sum of ring members, MEMB
c         changed subroutines: IsItCyclic and DefLoc 
c         (new lcs for atoms belonging to two planar rings)
c         X-H distances updated according to allen et al 2010
c         changed subroutine: ExtendH
c         other small cosmetic changes to make lsdb work properly
c         any probles listed out in file problems2010  
c 091009: Added exceptions in DetExcept for:
c         Nitrogen: N-c3conf
c         Oxygen: aromaturea, aromatketone, aromatester, aromatamide, urea
c 081114: Added exceptions in DetExcept for n2-C-o, Cnoc, C-n3, C-n2, c-N=o,
c-        N-oh, n-O, n-O2, c3pO, pO4conf, Po4conf
c 080808: Subroutines have been sorted in .f file
c 120206: O-H distance in phosphate changed to 1.015 
c 101006: Misc fix in printout when verbosity level >= 3
c 100606: Fix in reading HFIX/AFIX instructions from SHELX 
c         Misc fixes
c 092906: Misc fix in ReadInp (renaming of atoms from SHELX was not working)
c         Misc fixes to conform g77 standard 
c         Modified the way atoms from SHELX are renamed, now use chan choose 
c         between straighforward atom counting and counting of atoms in each
c           atom type 
c 081406: O-H distance in water changed to 0.9584
c         Misc fix for x-o-H in DetEx
c 080306: Fixed Ti --> Tl
c 073106: Fixed bug with ring search in several separate fragments
c         Added option to change default covalent radii of atoms
c             RADII Ca 0.04  O 0.45  H 0.59  N 0.49
c         Default covalent radii reset to original LSDB 
c         Added verbosity level
c             VERBOSE ilevel
c          ilevel = 0 : minimum info printout (really short)
c                 = 1 : only necessary info is printed [DEFAULT]
c                 = 2 : more info printed
c                 = 3 : even more info printed etc...
c         Added option of using line continuation character '\', but 
c           single instructuion line still can not exceed 255 characters  
c         A lot of misc fixes  
c 072806: changes by PMD in DetEX & ExtendH from PMD 072006 version of  
c         source code included
c 072806: Fixed serious bug in ring search...
c 072706: Planarity tolerance split into ring (tol2) and group (tol3)
c           TOLS bond 0.4 ring 0.1 group 0.1
c         Misc. fix in ExtendH
c         Removed interactive input, i.e. only input via input file
c           lsdb.inp is now supported
c 072606: Misc fix in ExtendH
c 072606: Added exception in DetEx for c2-C=o (ketone), c-Nc1-c1 N-planar
c 072506: Added creation of 'const.bond' file with CON instructions
c          for H-atoms when EXTEND_H option is activated
c 072106: Misc fix of reading scale factor from XD parameter file
c         Fixed issues with reading of XD files
c         Fixed issue with too small angle in X-Y-Z conformation 
c         Added automatic generation of 'RESET BOND' instructions when
c            extending of H-atom is requested 
c         Fixed issue with lmx(i) when reading XD files  
c         A lot of misc fixes 
c 072006: Added new subroutine to generate dummy atoms with checking all
c           previous dummy atoms in the list 
c         Subroutine DummyXYZ replaced with MidPoint
c         FIXED HORRIBLE BUG IN EQUIVALENCY OF KAPPA DETERMINATION !!!
c         Disorder (if any) is correctly applied to Pv's and Plm's of 
c           atoms, so correct sum of Pv (after DB) should be obtained
c         Added input option for selection of keys in key table
c           KEYGEN (*)XYZUIJ (*)all (*)ordered (*)MULT (*)all (*)ordered
c           XYZUIJ - atomic positions and thermal parameters 
c                   (Uij for heavy atoms and Uiso for hydrogens)
c               *all - for all atoms in the structure
c               *ordered - for ordered atoms ONLY
c           MULT - multipole parameters    
c               *all - for all atoms in the structure
c               *ordered - for ordered atoms ONLY
c         Misc fix in writing XD parameter file
c         A lot of misc fixes 
c 072006: Added exceptions in DetEx for cPo2(o-x) cPo(o-x)2 cPo3
c         cpO3 cpO2(o-x) cpO(o-x)2 x3-p-oH
c 072006: O-H in x-po3h & x-po3h2 extended to 1.015 Angstroms
c 071906: Misc fixes in handling of disorder disorder
c           DISORD  extendH  db
c           activate keywords if extending of disordered H-atoms to 
c           neutron distances and assignment of pseudoatom parameters 
c           to disordered atoms are desired (default is DO NOT for both)
c 071906: Added exceptions in DetEx for protonated aminophenyl
c 071806: Misc fixes of disorder in rings
c 071406: Damn! the ring search is still bad... now rewrote it
c           according to original papers of Corey & Wipke but it
c           seems to have problems with >2 fused ring systems
c         Fixed reading SHELX files with FVAR/HFIX/AFIX instructions 
c         Added calculation of Ueq
c         Added option for SHELX files to just create XD files bypassing
c           everything else (keeps OSF, WGHT etc)
c         Make sure that neighbors that are too close (i.e. 
c           < 0.5 Angstroms) and duplicate neighbors are not included 
c           (this often happens in disordered structures)
c         Misc fixes
c 070606: Fixed bug with renaming of SHELX atom names
c         Possibly (yet again!) fixed bug with ring search
c 042706: Added exceptions in DetEx for C & O atoms in carbamid ester
c 042006: Added exceptions in DetEx for n-=Cx=-n atoms
c 041706: Added exceptions in DetEx for n-=Cx=-c atoms
c 033006: Misc fix with 'iSimple' - ignore it when define kappa sets
c 021706: Misc fix in SameAtoms 
c         Fixed treatment of chemical equivalency of H-atoms
c         Fixed treatment of o-H H-atom of hemiacetal group in ExtendH
c 021506: Changes to ring perception algorithm by Marc Messerschmidt -
c           seems like it's working now :)
c         S-H extended to 1.338 Angstroms
c         Renaming atoms from SHELX is done by using keyword *RENAME 
c           under INPUT line
c         Determination of CHEMCON and KAPPA sets is based on SameAtoms
c           function and H-atoms are treated as other atoms
c         Added new option for using/not_using chemical/kappa constraints 
c           for H/non-H atoms
c              CHEMCON  *hydrogen  *other
c 021006: Added exceptions for H-atoms
c         A lot of additional exceptions in DetEx for C-atoms
c 020906: Misc fix of new algorithm of using sp-hybrid. in atomic environment
c         Fixed problem with 3m in non-planar sp2 group
c         Put in new exceptions (a lot of them!!!)
c 020806: Fix in special treatment of the ring atoms
c 020706: In SameAtoms add parameter which makes comparison much simpler
c         Added 'special' sorting of neighbors in FindAllNeigbors
c         New keyword 'special' in line 
c              RINGS [n1] [n2] (*)special 
c         determines if atoms in single planar rings are treated as 
c         usual atoms (Marc/Paulina) or a special case, when x-axis
c         is directed towards the center of the ring on other is directed
c         towards isn(1)
c 012006: Fixed small bugs in DetEx
c 011906: Fixed small bugs in DetEx
c 011706: Fixed small bug in DetEx w/carboxyl.. group introduced in 
c         previous version
c         For AMIDE group N just needs to be 3-coordinated
c         Fixed issues with systems that are not clear ETHERS/ANHYDR.etc
c 011306: Fixed issue with ESTER/ANHYDRIDE O-atoms in DetEx
c 101705: Misc improvements
c 092105: Added specification of NET CHARGES of fragments to lsdb.inp
c         Misc improvements
c 091405: Misc fix with 'isca'
c         Important fixes in 'DetExcept' and 'PrintGroup'
c 090705: Misc fix : asym = '1' --> 'NO'
c 090505: Misc fix for scaling of Pv's when NOT using DB (diff = default)
c 090105: Removed subroutine 'Save Nei' - not needed anymore
c         Fixed leftovers of 'SaveNei'
c 083105: Important fixes in ring perception - runs slightly faster now...
c         Databank file is opened only once in the very beginning
c         Added timings
c         Removed unused arrays of ~200MB !!!!!!!
c         Significant modifications to how the neighbors are handled :
c           Covalently-bonded neighbors (including sorting on distances 
c           and atom types)are precaluated for all atoms and stored in 
c             common block /neighbors/ 
c           'FindNei' needs to be called only when have to find 
c             non-covalently bonded neighbors  
c         SIGNIFICANT IMPROVEMENTS IN SPEED !!!  
c 083005: Fixes in 'DetExcept' etc...
c         Fixes in 'ReadLine' to make LSDB read properly the input 
c           file when compiled with INTEL Fortran
c         Got rid of logical variable 'ucov'
c         Added option of using connectivity matrix in 'Get Nei', so
c           that actual local system definition part runs MUCH FASTER!!
c         A lot of misc. fixes
c 082305: Implemented test for atoms in special positions
c         Now we can separate those from disordered atoms...
c         A lot of misc. fixes
c         Added specification of tolerances in input file
c         Added 'lsdb -h' option which prints out sample lsdb.inp file
c 082205: Yet another revision of RING PERCEPTION - seems to work fine!!!
c         A lot of misc. fixes
c         Removed command line arguments - created input file lsdb.inp
c 081905: Completely rewritten ring perception part - MUCH FASTER
c         minor changes on how to use Planarity tests etc.... 
c 081705: MAJOR CHANGES!!!! 
c         1. For crystal systems COMPLETE fragments are now generated
c            before entering search for rings
c         2. Neighbor handling subroutines now work on atomic sequence
c            numbers in common /atoms/
c         3. Things are much more simplified now...
c         4. If local system is defined with one of the sym.-gen. atoms
c            then at the end we rename it to a dummy atom...
c         5. 'DetExcept' is re-written
c         6. A lot of misc. fixes....
c         7. Create xd_lsdb.inp file if extending H-atoms
c 081105: automatic elongation of X-H bonds according to neutron values
c         misc fixes 
c 072805: important fix of dummy atoms counter for the case when
c         atoms which belong to one ring are NOT listed in a row!!!
c 072705: misc fix of symmetry label when reading databank entries
c         misc fix when reading databank file name 
c 030405: misc fix with connectivity when disordered atoms are present.
c         of course, it was Marc Messerschmidt who found the bug! 
c 022405: misc fix in ReadDB (default for lPGi was not correct)
c 022305: misc fixes (by default use ALL in sum-of-Pv rescaling)
c 022205: fixed bug with visualization of 3m (Z was wrong)
c 021805: added separate rescaling of several fragments - atoms
c           which belong to fragments can be specified in file
c           'frags.atoms' in the following format:
c            FRAG 1   1 -3  6 -13 18 -26
c            FRAG 2   4  5 14 -17 27 -36
c         default location of databank can be provided via environmental
c         variable LSDB_DB, i.e.
c           setenv LSDB_DB /crystal1/che9992/databank/organics-35.db
c           the above line can be also added to $HOME/.cshrc file
c         **** TODO: ADD THESE OPTIONS TO COMMAND LINE ****
c 021705: added test for planarity of 'bonded atomic group' of an atom
c         modified how the subroutine 'planarity' is called
c         added keyword 'PLANAR' to the databank
c         added support for MODEL name when reading XD parameter file
c 021605: misc fixes (ndum and reading rings info from command line)
c 121404: misc fixes to make all subroutines more portable
c         removed subroutine 'CleanChar'
c 120904: fixed ester group for O-atoms
c         slightly modified ReadDB when dealing with exceptions
c 102803: misc fixes
c 101303: misc fix in "DetExcept"
c 100603: subroutine for distinguish. except. groups "DetExcept" rewritten 
c         it seems to be working very well :)
c 093003: added exceptions for distinguishing similar atoms in the databank
c         and also for chemical equivalency. A lot of improvements!!!!
c 091503: fixed a bug with lattice types for SHELX
c 091203: added SIGMA (or ESD) of Pv to DataBank file and rescaling of
c         Pv's based on S. Price's method
c 070803: added optional EXCEPTION in DataBank format  (uh, ugly!!!)
c         (still needs some work!!!)  
c 070303: added scaling of Pv to any specified value
c         added two different method of scaling: multiplication and addition
c 022603: fixed reading PLATON file
c 022103: fixed bug with fractional<->cartesian coordinates of dummy atoms
c 021303: misc fix for linear conformation
c 021203: misc fix for dummy atoms in rings
c 020503: misc fix for X - ? - Y   group
c 013003: misc fix for 3m in sp3 atom
c 011703: smart decision on lmax for each atom when writing xd_lsdb.inp
c            - NOT good - it's not recognized by neither XD nor XDINTER
c 011603: symmetry codes '1' & '-1' replaced by 'NO' & 'C', respectively
c 011503: mix fixes in ReadDB, ReadChar, ReadLine, ReadFile 
c 011103: subroutine SYMOP is rewritten, or rather the "rotation" part 
c           of it. it should be much more rigorous now
c         fixed unit for SHELX file in ReadSHELX according to Louis's 
c           notes
c 011003: major revision - included symmetry elements when generating
c           neighbors of unique atoms
c         most subroutines require cartesian coordinates instead of 
c           references to unique atoms
c         NEEDS MUCH MORE DEBUGGING AND TESTING !!!!!!!!
c 010803: misc fixes by Piero & Louis
c         g77 compiles program just fine
c         misc cleanup of the code 
c 010703: misc fixes by Louis
c         fixed reading SHELX file
c         fixed dummy atom records in PLATON plot file 
c 010603: added support for dummy atom in PLATON plot of the local axes
c         misc fix for transformation matrices 'rm' and 'rm1'  
c         added subroutine 'DummyXYZ' - midpoint between two atoms 
c         misc cleanup of the code
c         misc fixes of PLATON instructions for local axes viz.
c 010303: misc fix in 'ReadPlaton' subroutine
c 121902: added SYMM keyword to databank file - additional test
c 121602: added PLATON plot with covalent radii
c 112702: misc fix (occupancies)
c 112602: when searching for chemical equivalence also make use of 
c           the symmetry itself (was not take in account before!)
c         misc improvement in ring comparison 
c         misc improvements in printout 
c 112002: misc improvement of the sequence of atoms in the ring
c 111402: yet another modification of ring perception subroutine,
c         hopefully this is the last one 
c 111202: ring perception subroutine is completely rewritted -
c          it is now sufficiently fast and accurate
c         connectivity table for ring perception is "reduced" 
c          by step-by-step elimination of open atoms
c 110602: many improvements
c 103102: added option to use Coppens Lab pseudoatom databank !!!!
c 092602: misc fix in 2+2 case, though I'm not sure it works properly now
c 061402: can read PLATON ANGSTROM format
c 021202: can read AIM2TAB output file
c 021102: misc fix in reading xd param file
c 012302: 'cyl.' --> 'cyl'
c         no key integers are written to chem.-equiv. atoms
c 012202: major modification of 'ReadShelx' subroutine - 
c         now it only asks for input file name and can identify
c         all SHELX keywords
c         added 'ReadXD' subroutine and unified output from
c          'ReadXD' and 'ReadShelx'
c         'ReadXD' can read 1 & 2 versions of xd parameter file
c         in case of reading XD files program will copy all other
c         keywords from old XD master file to xd_lsdb.mas
c 012102: changes in subroutine 'FindRings' - list of atoms
c         in the ring is now sorted on the connectivity
c         major fixes in inei=4 part
c         misc fixes everywhere
c         major bug in determination of kappa set - FIXED
c 011802: modified ring search - uses connectivity matrix
c         planarity test is now in separate subroutine  
c         misc fixes
c 011702: subroutine for H merged with subroutine for C
c         subroutine for O merged with subroutine for C
c         new subroutine is called "DefLoc"
c 011602: also figures out the kappa sets and chemical constraints
c
c 040520: Fixing/changing PLM (3,-3) & (4,-3) sign for mopro par file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c TODO: 
c -----
c   Treatment of disorder with 0.5/0.5 occupations
c   For periodic 3D networks stop growth of fragments in CompleFrags 
c      when sym.-generated atoms are different by just a translation
c   Add geometrical criteria for some higher symmetries, like mm2 (?) 
c     for cases with 4 and more neighbours
c   Group planarity tests for 4-coordinated atoms
c   ---------------------------------------------------------------
c   FORTRAN-90 : specification of dimension of arrays in input file
c   Add CIP R/S tests (see Platon)
c   XD WILL NOT READ CORRECTLY THE SYMMETRY CODES 2,222,3,4 etc
c   Add support for higher local symmetries (?)
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compilation: 
c ------------
c
c Lahey/Fujitsu Fortran95:
c   lf95 lsdb.f -o lsdb --nchk --trace -O --tpp --staticlink
c   lf95 lsdb.f -o lsdb  --chk --trace -g --staticlink
c   lf95 lsdb.f -o lsdb  --chk --trace -g --chkglobal --staticlink
c
c GNU Fortran77
c   g77 lsdb.f -o lsdb -fno-backslash -O2 -static
c   g77 lsdb.f -o lsdb -fno-backslash -Wall -Wsurprising -static
c
c INTEL Fortran
c   ifort  lsdb.f -o lsdb -O3 -tpp6
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prints out program header  
      subroutine Header(io)
      character line*71
      do i=1,len(line)
       line(i:i)='='
      enddo
      write(io,800)
800   format(/1X,71('=')/3x,'PROGRAM  L S D B  -  AUTOMATIC SETUP OF THE
     * LOCAL COORDINATE SYSTEM',/,1x,71('=')/,1x,
     *'Program determines the local coordinate system for each atom base
     *d on',/,16x,'connectivity and geometrical parameters',/,
     *' In addition it provides an interface to COPPENS LAB pseudoatom d
     *atabank',/,1x,
     * 'Sample LSDB input file (lsdb.inp) can be printed out with :',
     * ' lsdb -h',/,1x,71('-'),/,
     *' Unauthorized use of this program is strictly prohibited and subj
     *ect for',/,
     *' prosecution under the penalty of laws of the Russian Federation 
     *and USA',/,
     *1x,71('=')/
     *' Written by Dr. A. Volkov w/contrib. from M. Messerschmidt &',/,
     *' K.N. Jarzembska, P.M. Dominiak, P. Kumar',/,
     *' Copyright (c) 2006, A. Volkov',12x,
     *'Last modified on October 2018',/,1x,71('='))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prints out sample lsdb.inp file on the screen
      subroutine PrintHelp(io)
99    format(100a)
      write(io,99) ' '
      write(io,99) '!'
      write(io,99) '! Example of LSDB input file (lsdb.inp) : '
      write(io,99) '!'
      write(io,99) 'FORMAT  *shelx xd aim2tab platon '
      write(io,99) '! for SHELX :'
      write(io,99) 'INPUT   l-dopa.shelx   rename 0 toxd '
      write(io,99) '! for XD : '
      write(io,99) '!INPUT   xd.mas xd.inp '
      write(io,99) '! etc : '
      write(io,99) '!INPUT file.inp '
      write(io,99) 
     * '! printout level (higher number gives more printout) :'
      write(io,99) 'VERBOSE  4'
      write(io,99) 
     * '! tolerances for bonds / planar rings / planar groups'
      write(io,99) 'TOLS bond 0.4 ring 0.1 group 0.1'
      write(io,99) '! modify covalent radii'
      write(io,99) '!RADII  H 0.1  He 0.2  Li 0.3 ....'
      write(io,99) 
     *'! number of cells in each direction for complete fragment search'
      write(io,99) '! comment it out or set to -1 to turn off search'
      write(io,99) 'NCELLS  1'
      write(io,99) '! how to handle disordered part'
      write(io,99) 'DISORD  extendH  db'
      write(io,99) '! generation of keys'
      write(io,99) 'KEYGEN  XYZUIJ all ordered  MULT all ordered'
c
      write(io,99) 
     * '! use simple model for determination of chemical equivalency'
      write(io,99) 'SIMPLE'
      write(io,99) 
     * '! use/do_not_use chemical/kappa contraints for H/non-H atoms'
      write(io,99) 'CHEMCON  *hydrogen   *other'
      write(io,99) '! extend H-atoms to neutron/standard distance '
      write(io,99) 'Extend_H'
      write(io,99) '! search for 3-8 membered rings'
      write(io,99) 'RINGS  3  8  *special'
      write(io,99) 
     * '! assign pseudoatom parameters from COPPENS LAB databank'
      write(io,99) 'DATABANK /crystal1/che9992/proteins/db43ed081406.db'
      write(io,99) '! select method for scaling of Pv''s '
      write(io,99) 'SCALE *sigma diff fract'
      write(io,99) '! number of fragments '
      write(io,99) 'NFRAG 1'
      write(io,99) '! fragment specifications'
      write(io,99) '!FRAG 1 1 -25'
      write(io,99) '!FRAG 2 26 -50'
      write(io,99) '! net charges of fragments'
      write(io,99) '!CHFRAG 0. 0.'
      write(io,99) ' '
      return
      end     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine Init
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
      common /mpoles/ Pv(nat), sPv(nat), Plm(nat,0:4,-4:4), iar(nat,11)
      common /kappas/ nkap,rkappa(nat,2),nkapsf(nat)
      common /kappaset/ itkset, kset(nat)      
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      character amodel*20
      common /amodel/ amodel
      logical ldb
      common /ldb/ ldb, ntot, nfound, nmissed
      common /tols/ tol1, tol2, tol3
      character group(nat)*40
      common /except/ iexp(nat), group
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      logical lring
      character fname(2)*60
      common /args/ fname, lring, isca, sumofPV
      common /forH/ isp(nat), ihp(nat), ihnei(nat)
      common /shelx1/ fvar(200), wght(6), ifv
      logical led, ldd
      common /disorder/ led, ldd
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      logical lkeys
      common /lkeys/ lkeys(7)
      logical lrs
      common /renameSHELXatoms/ lrs, irs
      logical lsx
      common /lsx/ lsx
      logical lfH
      common /lfH/ lfH
      common /ifra/ ifra
      common /compfrag/ ncells
      logical lhyd, lheav
      common /lchem/ lhyd, lheav
      logical lspe
      common /ring1/ nrmin, nrmax, lspe
      common /ftype/ iform
      character abank*60
      common /abank/ abank
      common /simple/ iSimple
      common /verbosity/ iverbose

c.....Tolerances:
      tol1 = 0.2D00  ! tolerance added to the sum of covalent radii
      tol2 = 0.1D00  ! ring sigma plane ( below - planar, above - not )
      tol3 = 0.1D00  ! group sigma plane ( below - planar, above - not )

      iverbose = 1
      ncells = -1   
      ifra = 1
      nrmin = 0
      nrmax = 0
      iform = 0
      iSimple = 0
      sumofPV = 0.00D00

      abank = ' '

      lspe = .false.
      lhyd = .true.
      lheav = .true.
      ldb = .false.
      lfH = .false.
      lsx = .false.
      lrs = .false.     ! do not rename SHELX atoms
      irs = 0 
      lring = .false.
      led = .false.     ! exclude disorder from extend_H processing
      ldd = .false.     ! exclude disorder from DB processing

      amodel = 'TEST'
      fname(1) = ' '
      fname(2) = ' '
      isca = 0
      nkap = 0
      ntot = 0
      nfound = 0
      nmissed = 0
      nor = 0
      itkset = 0
      ifv = 0
c
      call rzero( wght, 6 )
      wght(1) = -2.00D00
      wght(6) = 1.D0/3.D0
      call rzero( fvar, 200 )
      fvar(1) = 1.00D00
c      
      call rzero( xf, nat*3 )
      call rzero( xc, nat*3 )
      call rzero( u, nat*6 )
      call rzero( Pv, nat )
      call rzero( sPv, nat )
      call rzero( Plm, nat*45 )
      call rzero( rkappa, nat*2 )
      call rzero( rnei, nat*nbond )
c
      call izero(   isp, nat )
      call izero(   ihp, nat )
      call izero( ihnei, nat )
      call izero(  nnei, nat )
      call izero(  jnei, nat*nbond )
      call izero( jtnei, nat )
      call izero( janei, nat*nbond )
      call izero( jnnei, nat*nbond )
      call izero( jsnei, nat*nbond*nbond )
      call izero( kset, nat)
      call izero( itype, nat )
      call izero( iut, nat )
      call izero( isft, nat )
      call izero( iaf, nat )
      call izero( lmx, nat )
      call izero( iar, nat*11 )
      call izero( nkapsf, nat )
      call izero( isymeq, nat)
      call izero( nsymeq, nat)
      call izero( iexp, nat)
      call izero( noar, maxrings )
      call izero( nra, maxrings*maxmem )
      call izero( ipldum, maxrings )
c      call izero(,)
c
      do i=1,nat
         atom(i) = ' '
         asym(i) = ' '
        group(i) = ' '
          occ(i) = 1.00D00
      enddo
c
      do i=1,maxrings
        planar(i) = .false.
      enddo
c
      do i=1,7
        lkeys(i) = .false.
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c this is the start of main program
      Program LSDB
      include 'lsdb.inc'
      dimension xj(3)
      character atom(nat)*8, asym(nat)*5
      character abank*60, line*255, ans
      character dlab*6,  mas*60
      character fname(2)*60
      character aa1*10, atom1*10, atom2*10
      character dstring*30
      character*90, dimension(nat,1) :: MOPRO
      logical ldb, lring
      common /mopro/ MOPRO
      character*100, dimension(nat,1)::RES
      common /cons/ RES
      character*17, dimension(nat,1):: RESI1
      character*10, dimension(nat,1):: RESI2
      logical resi11
      common /resi1/ RESI1, RESI2, resi11
      character*90, dimension(nat,1) :: MCON
      common /mcon/ MCON
      common /args/ fname, lring, isca, sumofPV
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /files/ imas, nout, ixp, idb, icon
      common /ftype/ iform
      common /ifra/ ifra
      common /iplat/ iplat
      common /ldb/ ldb, ntot, nfound, nmissed
      common /kappas/ nkap,rkappa(nat,2),nkapsf(nat)
      common /mpoles/ Pv(nat), sPv(nat), Plm(nat,0:4,-4:4), iar(nat,11)
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /abank/ abank
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      logical lspe
      common /ring1/ nrmin, nrmax, lspe
      common /at_types/ ntyp, iatype(nat)
      common /at_uij/ u(nat,6)
      common /covrad/ covrad(54)
      logical lfH
      common /lfH/ lfH
      character symbat(92)*2
      common /atsym/ symbat
      common /xdold/ mas
      common /nsft/ nsft
      character amodel*20
      common /amodel/ amodel
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      common /compfrag/ ncells
      logical lsx
      common /lsx/ lsx
      common /shelx1/ fvar(200), wght(6), ifv
      common /xhbonds/ irb, icb, ires
      real*4  yt, tm(2), tmt(2)
      data zero /0.00D00/
c---
c---  printout of the sample lsdb.inp file
c---
      if( iargc().ne.0) then
        amodel = ' '
        call Header(6)
        call getarg(1,amodel)
        if(amodel(1:2).eq.'-h') call PrintHelp(6)
        stop
      endif  
c---
c---  IO unit numbers
c---
c.....new master file (xd_lsdb.mas)
      imas = 58  
c.....new mopro file (mopro_lsdb.par)
      ipro = 444
c.....new mopro constraint file (res_mopro.inp)
      ires = 445 
c.....new mopro restraint file (con_mopro.inp)
      icon = 443 
c.....output/log file (lsdb.out)
      nout = 35
c.....PLATON file in ANGSTROM format file with local axes (lsdb.spf)
      iplat = 46      
c.....XD parameter file (xd_lsdb.inp)
      ixp = 48
c.....DATABANK file
      idb = 47
c.....misc SHELX-type file with H-positions moved to standard distances
      ifh = 49      
c.....file 'reset.bond' with RESET BOND instructions
      irb = 43
c.....file 'const.bond' with CONSTRAINTS instructions
      icb = 42
c.....also unit=52 - PLATON check file in ANGSTROM format in subroutine LocalDef
c                    (lsdb.ang)
c.....also unit=66 - PLATON check file in ANGSTROM format for reduced structure 
c                    in 'Connect4Rings' (reduced.ang)      
c.....also unit=57 - SHELX input file in 'ReadSHELX'
c
c--- initialize some arrays
      call Init
c---  
c---  print header
c---  
      io=6
      call Header(io)
c---                
      open(imas,file='xd_lsdb.mas',form='formatted',status='unknown')
c---                
      open(ipro,file='mopro_lsdb.par',form='formatted',status='unknown')
c---  header for mopro.par
        write(ipro,99) ('!',i=1,80)
        write(ipro,'(A)')'!               CREATED BY PROGRAM << LSDB >>'
        write(ipro,'(A)')'!                modified by K.N.Jarzembska'
        write(ipro,'(A)')'!                       Nancy, 2012'       
        write(ipro,99) ('!',i=1,80)
        write(ipro,99)'!  '   
c---                
      open(icon,file='con_mopro.inp',form='formatted',status='unknown')
c---  header for mopro.par
        write(icon,'(A)')'!'
        write(icon,'(A)')'! local symmetry  '
        write(icon,'(A)')'!'
c---      
      open(nout,file='lsdb.out',form='formatted',status='unknown')     
      call Header(nout)
c---  read and decode the arguments
      call ReadInp
c---  read structure files
      if(iform.eq.1) call ReadShelx( fname(1) )        
      if(iform.eq.2) call ReadXD( fname(1) )
      if(iform.eq.3) call ReadAim2Tab( fname(1) )
      if(iform.eq.5) call ReadPlaton( fname(1) )
c---  figure out the number of SFAC 
      nsft = 1
      do i=1,ia
        if( isft(i).gt.nsft ) nsft = isft(i)
      enddo   
c---  Platon plot file with local axes
      open(iplat,file='lsdb.spf',form='formatted',status='unknown')
      write(iplat,'(a)') 
     * 'TITL X-green, Y-yellow, Z-brown, Dummy(D)-orange'
      write(iplat,'(a)') 'ANGSTROM'
c     try to open databank
      if( .not. ldb ) goto 58
55    open(idb,file=abank,form='formatted',status='old',err=56)
      goto 58
56    write(*,'(/3a)') ' Error opening databank file ',abank
57    abank=' '
      write(*,'(1x,a,$)') 'Enter the name of databank file : '
      read(*,99) abank
      if(abank.eq.' ') goto 57
      goto 55
58    continue       

c---- for fractional coordinates - generate complete molecule(s) 
c---- by applying symmetry operations
      if( ifra.eq.1 .and. ncells.gt.-1 ) then
        call GenCompleteFrags
      else
        natoms = ia
      endif  
c----
c---- RINGS !!!!
c----
      if( lring ) then
        yt = dtime(tm)
        tmt(1) = tm(1)
        tmt(2) = tm(2)
c       get 'reduced' connectivity matrix:
        call Connect4Rings
c       actually find rings:         
        call FindRings 
        yt = dtime(tm)
        tmt(1) = tmt(1) + tm(1)
        tmt(2) = tmt(2) + tm(2)
c        yt = tm(1) + tm(2)
        write(nout,408) tm(1), tm(2), yt
        write(   *,408) tm(1), tm(2), yt
      endif
408   format(/' Timing for ring search u,s,t (s) :', 3(1x,f12.1) )
c---  
      if( lfh ) then
        open(irb,file='reset.bond',form='formatted',status='unknown')
        open(icb,file='const.bond',form='formatted',status='unknown')
      open(ires,file='res_mopro.inp',form='formatted',status='unknown')
        write(445, '(A)')'!' 
        write(445, '(A)')'! X-H distances' 
        write(445, '(A)')'!' 
      endif        
c---
c---  determine and store neighbors of each atom
c---
      call FindAllNeighbors
c---
c---  define local coordinate system      
c---
      call LocalDef
c---
      if( lfh ) then
           write(445,'(A)') '!'
           write(445,'(A)') '! isotropic H atoms ADPs'
           write(445,'(A)') '!'
           do i=1, nat
           if(RES(i,1)(1:6).eq.'URATIO') then
           write(445,'(A)') RES(i,1)
           endif
           enddo
        close(irb)
        close(icb)
        write(445,'(A)') '!' 
        write(445,'(A)') '! end'
        write(445,'(A)') '!'          
        close(ires)
      endif  
c---
c---  create new xd parameter file
c---
      if( ldb .or. lfH .or. lsx ) then
        open(ixp,file='xd_lsdb.inp',form='formatted',status='unknown')   
        write(nout,'(/80a)') ('-',i=1,80)
        write(nout,388)
        write(   *,388)
388     format(/' Creating XD parameter file...')
c---
        if( lsx ) then
c?????
cc         make 'lmx' of all atoms = 0
c          call izero( lmx, nat )
          goto 386
        endif  
c---
        if( lfH ) open(ifh,file='fixedH.ins',form='formatted')
c---    statistics
        if( ldb ) then
          write(   *,387) ntot, nfound, nmissed
          write(nout,387) ntot, nfound, nmissed
        endif   
387     format(/' ATOMS SEARCHED in DB, FOUND in DB, ',
     *          'NOT FOUND in DB = ',3i8)
c---    RESCALING POPULATIONS
        if( ldb .and. isca.ne.0 ) call rescalePv
c---    header for xd.inp
386     continue
        write(ixp,99) ('!',i=1,80)
        write(ixp,99)'!                   CREATED BY PROGRAM << LSDB >>'
        write(ixp,99) ('!',i=1,80)
        write(ixp,99)'XDPARFILE VERSION  2'
        write(ixp,99)'  ',amodel(:ils(amodel)),'   MODEL  4  2  1  0'
        write(ixp,99) 
     * 'LIMITS nat 300 ntx 31 lmx 4 nzz 30 nto 0 nsc 20 ntb 20 nov 2500'
        write(ixp,95) ia, nkap, nsft, ndum
        write(ixp,96) (zero,i=1,14)
        do i=1,ndum
          write(ixp,'(3f13.6)') (xdum(i,j),j=1,3)
        enddo
ccccc to MOPRO
      write(444,'(A,I3)') 'DUMMY ',ndum
      do i=1,ndum
        write(dstring,'(I5)') i
        write(444,'(A,3F12.6,A,A)') 'DUMY ',(xdum(i,k),k=1,3),
     *'  D'//adjustl(trim(dstring)),' 0'
      end do
      write(444,'(A)') ''
c---  kappa sets for MoPro
      write(444,'(A,I3)') 'KAPPA ',nkap
      do i=1,nkap
        write(444,'(2f10.6)') (rkappa(i,k),k=1,2)
      enddo 
      write(444,'(A)') ''
      write(444,'(A,I4)') 'ATOMS ',ia
      write(444,'(A)') ''
cccccccccc
95      format('USAGE ',i4,' 2 4 ',i4,' 0 1 0 0 ',i4,' 0 0 0 0 ',i4)       
96      format(8f10.6)
c--- atom entries in xd parameter file
        do k=1,ia
          write(ixp,33) atom(k),(iar(k,j),j=1,11),(xf(k,j),j=1,3),occ(k)
          write(ixp,'(6f10.6)') (u(k,j),j=1,6)
c          lmax=iar(k,9)
c          if(lmax.ge.2) write(ixp,34) 
          write(ixp,34) 
     *      Pv(k), Plm(k,0,0), Plm(k,1,1),Plm(k,1,-1), Plm(k,1,0),
     *      Plm(k,2,0),Plm(k,2,1),Plm(k,2,-1),Plm(k,2,2),Plm(k,2,-2), 
     *      Plm(k,3,0),Plm(k,3,1),Plm(k,3,-1),Plm(k,3,2),Plm(k,3,-2),
     *      Plm(k,3,3),Plm(k,3,-3),Plm(k,4,0),Plm(k,4,1),Plm(k,4,-1),
     *      Plm(k,4,2),Plm(k,4,-2),Plm(k,4,3),Plm(k,4,-3),Plm(k,4,4),
     *      Plm(k,4,-4) 
          if( lfH ) write(ifh,32) atom(k), isft(k), (xf(k,j),j=1,3), 
     *                            occ(k), u(k,1)
ccccomment
          atom1=atom(k)
          call erbras(atom1)
          atom2=atom(k)
          call ernum(atom2)
          if(resi11) then
          write(444,'(A,I3,1X,A,3F12.6,F6.3,1X,A,A)') 'ATOM ',k, 
     *    RESI1(k,1),(xf(k,j),j=1,3), occ(k), ' 1 ',atom2 
          else
          write(444,'(A,I3,1X,A,A,3F12.6,F6.3,1X,A,A)') 'ATOM ',k, 
     *    atom1,'1   MOL',(xf(k,j),j=1,3), occ(k), ' 1 ',atom2 
          endif     
          write(444,99) MOPRO(k,1)        
          aa1=atom(k)
          call ernum(aa1)
          if( aa1.eq.'H')then
          write(444,'(A,f10.6)') 'UISO   ',(u(k,1))
          else
          write(444,'(A,6f10.6)') 'UANI   ',(u(k,j),j=1,6)
          endif
          if( aa1.eq.'H')then
          write(444,'(2f8.5,2x,8(f5.3,2x))') 
     *      Pv(k)/occ(k), Plm(k,0,0)/occ(k), Plm(k,1,1)/occ(k),
     *      Plm(k,1,-1)/occ(k), Plm(k,1,0)/occ(k),
     *      Plm(k,2,0)/occ(k),Plm(k,2,1)/occ(k),
     *      Plm(k,2,-1)/occ(k),Plm(k,2,2)/occ(k),Plm(k,2,-2)/occ(k)  
          write(444,'(A)') ''      
          else
          write(444,'(2f8.5,2x,8(f5.3,2x))') 
     *      Pv(k)/occ(k), Plm(k,0,0)/occ(k), Plm(k,1,1)/occ(k),
     *      Plm(k,1,-1)/occ(k), Plm(k,1,0)/occ(k),
     *      Plm(k,2,0)/occ(k),Plm(k,2,1)/occ(k),
     *      Plm(k,2,-1)/occ(k),Plm(k,2,2)/occ(k),Plm(k,2,-2)/occ(k)
          write(444,'(2x,7(f5.3,2x))')         
     *      Plm(k,3,0)/occ(k),Plm(k,3,1)/occ(k),Plm(k,3,-1)/occ(k),
     *      Plm(k,3,2)/occ(k),Plm(k,3,-2)/occ(k),
cPK  PLM 3 -3 sign changed (multiplied with -1) Christain/PMD 03May2020
     *      Plm(k,3,3)/occ(k),(Plm(k,3,-3)/occ(k))*-1
          write(444,'(2x,9(f5.3,2x))')
     *      Plm(k,4,0)/occ(k),Plm(k,4,1)/occ(k),Plm(k,4,-1)/occ(k),
     *      Plm(k,4,2)/occ(k),Plm(k,4,-2)/occ(k),Plm(k,4,3)/occ(k),
cPK  PLM 4 -3 sign changed (multiplied with -1) Christain/PMD 03May2020
     *      (Plm(k,4,-3)/occ(k))*-1,Plm(k,4,4)/occ(k),
     *      Plm(k,4,-4)/occ(k) 
          write(444,'(A)') ''  
          endif
          enddo 
cccomment
           write(443,'(A)') '!'
           write(443,'(A)') '! similarity (switched-off)'
           write(443,'(A)') '!'
           do i=1, nat
           if(MCON(i,1)(1:8).eq.'! CONPVM') then
           write(443,'(A)') MCON(i,1)
           endif
           enddo   
           write(443,'(A)') '!'
           write(443,'(A)') '! end'
           write(443,'(A)') '!'  
           close(icon) 
cccomment 
           write(444,'(A)')''
           write(444,'(A)')'ANHAR    0    3         !  #atoms #order '
           write(444,'(A)')'! 1         ! atom number '
           write(444,'(A)')'! 0. 0. 0. 0. 0. 0. 0. 0. '
           write(444,'(A)')'! 0. 0. '
           write(444,'(A)')''
           write(444,'(A)')'STOP '
33      format(a8,2i3,3i5,2x,6i5,3f9.5,f8.4)
34      format(10f8.4)
32      format(a,1x,i4,2x,3f12.6,2x,f10.5,1x,f8.4)
c---    kappa sets
        do i=1,nkap
          write(ixp,'(i3,6f10.6)') nkapsf(i), (rkappa(i,k),k=1,2), 
     *    rkappa(i,2),rkappa(i,2),rkappa(i,2),rkappa(i,2)
        enddo 
        write(ixp,'(7e11.4)') (fvar(i),i=2,9)
        write(ixp,'(e13.6)') fvar(1)
        close(ixp)
        if( lfH ) close(ifh)
      endif  
c---- add dummy atoms to PLATON file (cartesian)
      do i=1,ndum
        if(ifra.eq.1) then
          call rabs(0,xdum(i,1),xdum(i,2),xdum(i,3),xj(1),xj(2),xj(3))
          write(iplat,35) i, xj    
        else
          write(iplat,35) i, (xdum(i,k),k=1,3)    
        endif  
      enddo      
35    format('ATOM D(',i1,')',2x,3f14.4)
c---- add plot instructions to PLATON file     
      if(ndum.eq.0) write(iplat,400) 
     *  ( symbat(iatype(k)) ,covrad(iatype(k)), k=1,ntyp ), 
     *   'X',0.01,'Y',0.01,'Z',0.01
      if(ndum.ge.1) write(iplat,400)
     *  ( symbat(iatype(k)) ,covrad(iatype(k)), k=1,ntyp ), 
     *   'X',0.01,'Y',0.01,'Z',0.01,'D',0.01
      write(iplat,401)
      if(ndum.ge.1) write(iplat,402)
400   format('JOIN RADII ',100(1x,a2,f5.2))      
401   format(
     *'SOLID COLOR',/,
     *'COLOR TYPE X GREEN Y YELLOW Z BROWN D ORANGE',/,
     *'RADII ATOMS ALL 0.06',/,
     *'RADII ATOMS X 0.02',/,
     *'RADII ATOMS Y 0.02',/,
     *'RADII ATOMS Z 0.02',/,
     *'RADII BONDS ALL 0.01 1',/,
     *'RADII BONDS TO X 0.01 1',/,
     *'RADII BONDS TO Y 0.01 1',/,
     *'RADII BONDS TO Z 0.01 1',/,
     *'RADII BONDS TO D 0.00 0',/,
     *'DETACH X Y',/,
     *'DETACH X Z',/,
     *'DETACH Y Z' ,/,
     *'DETACH X X',/,
     *'DETACH Y Y',/,
     *'DETACH Z Z' ,/,
     *'LABEL ON',/,
     *'LABEL H',/,
     *'UNLABEL X Y Z')
402   format(
     *'RADII ATOMS D 0.02',/,
     *'DETACH D X' ,/,
     *'DETACH D Y' ,/,
     *'DETACH D Z' ,/,
     *'LABEL D')
      if(ndum.ge.1) then
        do k=1,ntyp
          write(iplat,99) 'DETACH D ',symbat(iatype(k))
        enddo
      endif  
      write(iplat,99) 'PLOT'     
      close(iplat)      
c---
      write(nout,'(/80a)') ('-',i=1,80)
      write(nout,389)
      write(   *,389)
      if( ldb .or. lfH .or. lsx ) then
        write(   *,390)
        write(nout,390)
      endif
      if( lfH ) then
        write(   *,393)
        write(nout,393)
      endif
      write(   *,391)
      write(nout,391)
389   format(/' XD instructions are in file            : xd_lsdb.mas')
390   format( ' XD parameters are in file              : xd_lsdb.inp')
393   format(' Atomic coordinates w/extended H-atoms  : fixedH.ins',/,
     *       ' RESET BOND instructions for H-atoms    : reset.bond',/,
     *       ' CONSTRAINT instructions for H-atoms    : const.bond')
392   format(' Summary of DB atom types : grep ''DB comment'' lsdb.out')
394   format(' Summary of X-H distances : grep -i ''x-h'' lsdb.out' )
395   format(' Summary of exceptions    : grep -i ''Exception',
     *       ' identified'' lsdb.out' )
391   format(     
     *  ' Cartesian coordinates in PLATON format : lsdb.ang',/,
     *  ' Reduced structure in PLATON format     : reduced.ang &',
     *                                             ' reduced_1.ang',/,
     *  ' PLATON plot file with local axes       : lsdb.spf',/,
     *  ' Program log file                       : lsdb.out',//,
     *  ' Local axes can be visualized in PLATON using command : ',
     *'platon -p lsdb.spf')
      if( ldb ) then
        write(   *,392)
        write(nout,392)
      endif                   
      if( lfH ) then
        write(   *,394)
        write(nout,394)
      endif  
      write(   *,395)
      write(nout,395)
      close(imas)
      if( ldb ) close(idb)
      
c     timing      
      yt = dtime(tm)
      tmt(1) = tm(1) + tm(1)
      tmt(2) = tm(2) + tm(2)      
      write(nout,409) tmt(1), tmt(2), tmt(1)+tmt(2)
      write(   *,409) tmt(1), tmt(2), tmt(1)+tmt(2)
409   format(/' Total program elapsed times u,s,t (s) :', 3(1x,f12.1) /)
      
c
      write(nout,'(80a/)') ('-',i=1,80)
      close(nout)
c---
98    format(/a,$)
99    format(100a)  
      stop
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c reads input file
c
      subroutine ReadInp
      include 'lsdb.inc'
c
      character line*255
      common /ftype/ iform
      common /ifra/ ifra
      common /files/ imas, nout, ixp, idb, icon
      logical lring
      character fname(2)*60
      common /args/ fname, lring, isca, sumofPV
      character mas*60
      common /xdold/ mas
      logical lfH
      common /lfH/ lfH
      logical ldb
      common /ldb/ ldb, ntot, nfound, nmissed
      character abank*60
      common /abank/ abank
      logical lspe
      common /ring1/ nrmin, nrmax, lspe
      logical lrs
      common /renameSHELXatoms/ lrs, irs
      logical lsx
      common /lsx/ lsx
      logical led, ldd
      common /disorder/ led, ldd
      logical lkeys
      common /lkeys/ lkeys(7)
      common /compfrag/ ncells
      common /tols/ tol1, tol2, tol3
      common /simple/ iSimple
      logical lhyd, lheav
      common /lchem/ lhyd, lheav
      common /covrad/ covrad(54)
      character symbat(92)*2
      common /atsym/ symbat
      common /verbosity/ iverbose
      dimension rdum(20)
      character adum*2
      character ach(0:1)*24
      data ach /'straightforward counting',
     *          'counting atoms types    ' /  
c
      parameter( maxchar = 40 ) ! max number of character variables in the line
      character vchar(maxchar)*8
      integer istar(maxchar), ise(maxchar,2)
      common /line_extract/ nch,istar,vchar
c
90    format(100a)
c
      line = ' '
      mas = ' '
c      
      ii = 23
      open(ii,file='lsdb.inp',form='formatted',status='old',err=70) 
c.....printout      
31    format(/1x,71('-'))
      write(nout,31)
      write(nout,'(/a/)') ' LSDB options : '
c     when using input file the default for scaling is DIFFERENCE 
      isca = 1
10    call ReadInputLine(ii,line,*30)
        call ReadLine(line)
c        write(*,90) ' vchar(1) = [',vchar(1),']'
c       FORMAT
        if(vchar(1).eq.'FORMAT') then
          do i=2,nch
            if(vchar(i).eq.  'SHELX'.and.istar(i).eq.1) iform = 1
            if(vchar(i).eq.     'XD'.and.istar(i).eq.1) iform = 2
            if(vchar(i).eq.'AIM2TAB'.and.istar(i).eq.1) iform = 3
            if(vchar(i).eq. 'PLATON'.and.istar(i).eq.1) iform = 5
          enddo
          if(iform.gt.2) ifra = 0
        endif
c       INPUT
        if(vchar(1).eq.'INPUT') then
          do i=nch-1,nch
            if( vchar(i).eq.'RENAME' .and. istar(i).eq.1 ) lrs= .true.
            if( vchar(i).eq.  'TOXD' .and. istar(i).eq.1 ) lsx= .true.
          enddo  
c         how to rename SHELX atoms 
          if( lrs ) then
            call DNumLine(line,inum)
            write(*,*) ' inum = ',inum
            if( inum.ne.0 ) then
              call ReadChar(line,inum,i1,i2)
              read(line(i1-1:),*) irs
            endif  
          endif
          backspace(ii)
          call ReadInputLine(ii,line,*30)
          call DNumLine(line,inum)
          ilast = 2
          if( iform.eq.2 ) ilast=3
          do i=2,ilast
            call ReadChar(line,i,i1,i2)
            fname(i-1) = line(i1:i2)
          enddo
          if(iform.eq.2) then
            mas = fname(1)
            fname(1) = fname(2)    
          endif
        endif        
c       VERBOSITY
        if(vchar(1).eq.'VERBOSE') read(line,*,err=78) iverbose
c       TOLERANCE PARAMETERS
        if(vchar(1).eq.'TOLS') read(line,*,err=78) tol1, tol2, tol3
c       Atomic radii
        if(vchar(1).eq.'RADII') then
          call DNumLine(line,inum)
          if( inum.ne.nch-1 ) goto 78 
          read(line,*) (rdum(k),k=1,inum)
          do i=2,nch
            adum = ' '
            adum(1:2) = vchar(i)(1:2)
            call LoCase( adum(2:2) )
            do k=1,54
              if( adum(1:2).eq.symbat(k)(1:2) ) then
                write(nout,4669) symbat(k), rdum(i-1)
                covrad(k) = rdum(i-1) 
              endif
            enddo
          enddo
        endif
4669    format(' Changing covalent radius of ',a,' to ',f8.3,
     *         ' Angstroms')        
c       TRY TO FIND COMPLETE FRAGMENTS WITHIN -'ncells' < cells < +'ncells'
        if(vchar(1).eq.'NCELLS') read(line,*,err=78) ncells
c       disorder
        if(vchar(1).eq.'DISORD') then
          if( istar(2).eq.1 ) led = .true.
          if( istar(3).eq.1 ) ldd = .true.
        endif
c       keygen
        if(vchar(1).eq.'KEYGEN') then
          do i=2,nch
           if( istar(i).eq.1 ) lkeys(i-1) = .true.
          enddo 
        endif
c       EXTEND H-atoms
        if(vchar(1).eq.'EXTEND_H') lfH = .true.
c       CHEMICAL CONSTRAINTS
        if(vchar(1).eq.'CHEMCON') then
          if( istar(2).eq.0 ) lhyd = .false.
          if( istar(3).eq.0 ) lheav = .false.
        endif   
c       Use simple way of determining the "equivalence" of atoms
        if(vchar(1).eq.'SIMPLE') iSimple = 1
c       USE DATABANK      
        if(vchar(1).eq.'DATABANK') then
          ldb = .true.
          backspace(ii)
          call ReadInputLine(ii,line,*30)
          call DNumLine(line,inum)
          if(inum.le.1) then
            call getenv('LSDB_DB',abank)           
          else
            call ReadChar(line,2,i1,i2)
            abank = line(i1:i2)
          endif
        endif
c       RINGS
        if(vchar(1).eq.'RINGS') then
          lring = .true.
          read(line,*,err=78) nrmin, nrmax
          if(istar(2).eq.1) lspe=.true.
        endif
c       SCALING
        if(vchar(1).eq.'SCALE') then
          do i=2,nch
            if(vchar(i).eq.'SIGMA'.and.istar(i).eq.1) isca = 3
            if(vchar(i).eq. 'DIFF'.and.istar(i).eq.1) isca = 1
            if(vchar(i).eq. 'FRAC'.and.istar(i).eq.1) isca = 2
          enddo
        endif
      goto 10  
30    continue
      close(ii)      
c.....errors
      if( iform.eq.0 ) call Err('ReadInp','Unknown input file format')
      if( isca.eq.0  ) call Err('ReadInp','Unknown scaling method')
c
      write(nout,'(a,i2)') ' Verbosity level = ',iverbose
c
      if(iform.eq.1) write(nout,32) 'SHELX', fname(1)(:ils(fname(1)))
      if(iform.eq.2) write(nout,32) 'XD', mas(:ils(mas)),
     *   ' ',fname(1)(:ils(fname(1)))
      if(iform.eq.3) write(nout,32) 'AIM2TAB', fname(1)(:ils(fname(1)))
      if(iform.eq.5) write(nout,32) 'PLATON', fname(1)(:ils(fname(1)))
32    format(' Input file format = ',a,4x,'file(s) = ',4a)
c
      if( lrs ) then
        if( irs.ne.0 .and. irs.ne.1 ) irs = 1
        ij = ils( ach(irs) )
        write(nout,33) irs, ach(irs)(:ij)
      endif  
33    format(' Change SHELX atom names ( irs=',i1,' : ',a,' ) ')      
c
      if( lsx ) then
        write(nout,34)
        write(nout,31) 
        lfh = .false.
        ldb = .false.
        lring = .false.
        return 
      endif  
34    format(/
     * ' SIMPLE CONVERSION OF SHELX FILE TO XD FILES WAS REQUESTED',/,
     * ' ALL OTHER OPTIONS WILL BE IGNORED !!!!!!')
c
      if( ncells.eq.-1 ) then
        write(nout,90) ' Do not search for complete fragments'
      else
        ncells=abs(ncells)
        write(nout,35) -ncells, ncells
      endif
35    format(' Search for complete fragments within ',i3,
     * ' ... ',sp,i3,' unit cells')
      call PrtTol(nout)
c
      write(nout,90) ' Treatment of disorder :'
      if(     led) write(nout,90) '   apply Extend_H option'      
      if(.not.led) write(nout,90) '   do NOT apply Extend_H option'      
      if(     ldd) write(nout,90) '   apply DB'      
      if(.not.ldd) write(nout,90) '   do NOT apply DB'  
      if( ldb.and.ldd .and. .not.led ) write(nout,40) 
40    format(3x,'WARNING!!! Applying DB to disordered H-atoms',/,
     *                 14x,'w/o extending distances is CRAZY !!!')
c
      write(nout,90) ' Instructions for KEY table generation:'
      if( .not.lkeys(1) ) then 
        lkeys(2) = .false.
        lkeys(3) = .false.
      endif 
      if( .not.lkeys(4) ) then
        lkeys(5) = .false.
        lkeys(6) = .false.
        lkeys(7) = .false.
      endif 
      if( lkeys(2).and.lkeys(3) ) lkeys(2) = .false.
      if( lkeys(5).and.lkeys(6) ) lkeys(5) = .false.
      write(nout,44) lkeys
44    format(3x,'XYZUIJ=.',l1,'. all=.',l1,'. ordered=.',l1,'.',/,
     *       3x,'  MULT=.',l1,'. all=.',l1,'. ordered=.',l1,
     *          '. kappas=.',l1,'.')
c
      if( ldb ) write(nout,36) abank(:ils(abank))
36    format(' Read pseudoatom parameters from databank file :',/,3x,a)
c
      if( iSimple.ne.0 ) write(nout,37) 
37    format(' Use simple model for determination of ', 
     *       'chemical equivalency')
c
      if(      lhyd ) write(nout,38) 'U', 'H-atoms'
      if( .not.lhyd ) write(nout,38) 'Do NOT u', 'H-atoms'
      if(      lheav) write(nout,38) 'U', 'non-H atoms'
      if( .not.lheav) write(nout,38) 'Do NOT u', 'non-H atoms'
38    format(1x,a,'se chemical/kappa constraints for ',a)
c
      if(lring) write(nout,39) nrmin, nrmax
39    format(' Find rings with ',i3,' - ',i3,' atoms') 
c
      if(lspe) write(nout,41) 
     *'always m-symmetry : X - towards center, Y - neighbor in the ring'
      if(.not.lspe) write(nout,41) 
     *'as other atoms (only mark with iring=1)'
41    format(' Treatment of atoms that belong to 1 planar ring :',/3x,a)
c
      if(lfH) write(nout,42)
42    format(' Extend H-atoms to standard/neutron distances and create',
     * ' files : ',/,
     *  3x,'''reset.bond'' with RESET BOND instructions',/,
     *  3x,'''const.bond'' with CONSTRAINTS instructions')      
c
      if( ldb ) then
        if( isca.eq.1 ) write(nout,43) 'difference'
        if( isca.eq.2 ) write(nout,43) 'fraction'
        if( isca.eq.3 ) write(nout,43) 'sigmaPv'
      endif  
43    format(' Scale total population using ',a,' method')
c
      write(nout,31) 
      write(nout,'(/a)') ' Atomic covalent radii (Angstroms) : '
      write(nout,404) ( symbat(k), covrad(k), k=1,54 )
404   format( /7( 3x,a2,f5.2 ) )
      write(nout,31) 
c
      return
c.....
70    call Err('ReadInp','Unable to open LSDB input file lsdb.inp')
c.....
78    call Err('ReadInp','Error while reading file lsdb.inp')
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read one line from file 'ii'
c ignore lines commented with '!'
c make sure that continuation of line '\' is supported
      subroutine ReadInputLine(ii,line,*)
      character line*(*), dumline*255
90    format(100a)
      line = ' '
c      write(*,*) len(line)
c      stop
      dumline = ' '
      i1 = 1      
10    read(ii,90,err=71,end=70) dumline
       if( dumline(1:1).eq.'!' .or. dumline.eq.' ' ) goto 10
       ij = ils(dumline)
       if( dumline(ij:ij).ne.'\' ) then
         line(i1:) = dumline(:ij)
         return
       else
c        check total length of 'line'
         if( (i1+ij) .gt. len(line) ) goto 72
         line(i1:) = dumline(:ij)
         i1 = ils(line) + 2       
         goto 10
       endif
c
70    return 1
71    call Err('ReadInputLine','Error reading input file')
72    call Err('ReadInputLine','Line exceeds limit of 255 characters')
c
      end      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prints out tolerances
      subroutine PrtTol(io)
      common /tols/ tol1, tol2, tol3
      write(io,456) tol1, tol2, tol3
456   format(' Threshold / Tolerance parameters : ',/,
     * 3x,'bond criterion (tol1) : d(i-j) < covR(i) + covR(j) + ',f6.3,
     *       ' Angstroms',/,
     * 3x,'ring  planarity criterion (tol2) : sigpln < ',f7.3,
     *    ' Angstroms',/,
     * 3x,'group planarity criterion (tol3) : sigpln < ',f7.3,
     *    ' Angstroms')
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c READ XD MASTER AND PARAMETER FILES
c
      subroutine ReadXD(par)
      include 'lsdb.inc'
      parameter( maxchar = 40 ) ! max number of character variables in the line
c
      character mas*60, par*60, line*255
      character atom(nat)*8, asym(nat)*5, latmod*1
      dimension idum1(22)
      character vchar(maxchar)*8
      integer istar(maxchar)
      common /line_extract/ nch,istar,vchar
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
      common /files/ imasnew, nout, ixp, idb, icon
      common /xdold/ mas
      common /at_types/ ntyp, iatype(nat)
      common /ucell/ uc(6)
      character amodel*20
      common /amodel/ amodel
      common /shelx1/ fvar(200), wght(6), ifv
c
      amodel=' '
      imas=11
      ipar=12
c      mas=' '
c      par=' '
      line=' '
c
c691   write(*,'(a,$)') 
c     *' Enter XD master filename [xd.mas] : '
c      read(*,'(a)') mas
c      if(mas.eq.' ') mas='xd.mas'
      open(imas,file=mas,form='formatted',status='old',err=691)
c692   write(*,'(a,$)') 
c     *' Enter XD parameter filename [xd.inp] : '
c      read(*,'(a)') par
c      if(par.eq.' ') par='xd.inp'
      open(ipar,file=par,form='formatted',status='old',err=692)
      icent=1 ! centrosymmetric lattice
      latmod='P' ! primitive lattice
c
c     read xd master file
c
      write(nout,'(/2a)') ' Reading data from XD master file: ',mas
40    read(imas,'(a)',end=41,err=693) line
        if(line(1:5).eq.'CELL ') read(line(5:100),*,err=693) uc
        if(line(1:5).eq.'LATT ') then
          call ReadChar(line,2,i1,i2)
          if(line(i1:i1+1).eq.'A') icent=0
          call ReadChar(line,3,i1,i2)
          latmod=line(i1:i1+1)
          call UpCase(latmod)
        endif
      goto 40
41    continue        
      write(nout,'(/a,3f10.3,3f8.2)') ' Cell parameters: ',uc
      write(nout,'(/2a)') ' Lattice type : ',latmod
      if(icent.eq.1) write(nout,42) 'CentroSymmetric'
      if(icent.eq.0) write(nout,42) 'Acentric'
42    format(' Lattice is ',a)
c---  symmetry operations
      write(nout,'(/a)') ' Symmetry operations : '
      call GetSym(imas,nout)
      if(icent.eq.1) call AddCentSymm(nout)
      if(latmod.ne.'P') call NonP(nout,latmod)      
      write(nout,'(/a)') 
     * ' Final set of symmetry operations in matrix form:'   
      call PrSymOp(nout)
c---  transformation matrices
      call TranMat
      close(imas)
c
c     read xd parameter file
c
      write(nout,'(/2a)') ' Reading data from XD parameter file: ',par
      line=' '
129   format(/a,i4)  
100   format(100a)        
c--- skip the header of xd file:
99    call ReadXDLine(ipar,line)
      call ReadLine(line)     
c.....first figure out what the XD parameter file version is:
      ipv=1
C      write(*,'(3a)') 'line [',line,']'
C      write(*,'(3a)') '[',vchar(1),']'
      if(vchar(1).eq.'XDPARFIL') read(line,*) ipv
      write(nout,129) ' XD parameter file version = ',ipv
      if(ipv.eq.2) then
        call ReadXDLine(ipar,line)
        call ReadLine(line)
      endif
c.....MODEL
c      write(*,'(3a)') 'line [',line,']'
c      pause
      amodel = vchar(1)
      read(line,*) i1, i2, i3, i4
      write(nout,129)' Maximum level of multipole expansion = ',i1
      write(nout,100)' Model name = ',amodel
c--- read first record of xd parameter file: 
      line=' '
389   call ReadXDLine(ipar,line)
      call ReadLine(line)    
      if(ipv.eq.1) then
        read(line,*) idum1
      else
        read(line,*) (idum1(k),k=1,8)
        call ReadXDLine(ipar,line)
        call ReadLine(line)
        read(line,*) (idum1(k),k=9,22)
      endif  
      ia = idum1(9)   
      nz = idum1(12)  
      write(nout,'(a,i4)')' Number of atoms in XD parameter file = ',ia      
      write(nout,'(a,i4)')' Number kappa sets = ',nz
      if(ia.gt.nat) call Err('ReadXD','Too many atoms')
      write(nout,345) 
      write(   *,345) 
345   format(/1x,'NO.',3x,'ATOM',6x,'AT.NO.',3x,'ORDER_OF_U',5x,'iSFT')      
c      if(ia.gt.nat) then
c        write(*,400) ia,nat
c400     format(/'0Number of atoms in XD file (',i4,') is greater than '
c     */,' program limit (',i4,'). Please recompile the program...',//
c     * ,'0Program terminated...'/)
c        stop
c      endif
c--- skip dummy atoms & stats:
      do i=1,idum1(22)+2
        read(ipar,'(a1)') line
      enddo
c--- read atoms:
      do i=1,ia
567     call ReadXDLine(ipar,line)
c        write(*,'(/3a)') 'line=[',line(1:ils(line)),']'
        call ReadLine(line)
c        write(*,'( 3a)') 'line=[',line(1:ils(line)),']'
        atom(i)(1:8)=vchar(1)  
        read(line,*,err=567) (idum1(j),j=1,11), (xf(i,k),k=1,3), occ(i)
        iut(i)=idum1(6)   ! order of displacement tensor   
        isft(i)=idum1(7)  ! scattering factor table 
c       don't like because what if we want to refine multipoles...
c        lmx(i)=idum1(9)  
c        write(*,*) ' i, iut = ',i,iut(i)
c---    read U's
        read(ipar,*) (u(i,j),j=1,6)
        igo=0    
c       skip higher U's
        if(iut(i).eq.3) igo=2
        if(iut(i).eq.4) igo=5
        do j=1,igo
          read(ipar,'(a)') line
        enddo         
c.......get atomic number
        itype(i) = iAtNum( atom(i) )        
c---    default multipole expansion lmax=4
        lmx(i)=4
c---    for H-atoms - lmax=2       
        if(itype(i).eq.1) lmx(i)=2     
c
        write(   *,869) i, atom(i), itype(i), iut(i), isft(i)
        write(nout,869) i, atom(i), itype(i), iut(i), isft(i)
869     format(i4,3x,a,2x,i3,3x,i8,5x,i8)
c---    skip multipoles - DOESN'T ALWAYS WORK PROPERLY !!!!
        npop=idum1(9)*(idum1(9)+2)+1
        iskip=max(1,int(npop/10)+1)
c        iskip=max(2,int(npop/10)+1)
c        write(*,*) 'lpole=',idum1(9),'  npop=',npop,' iskip=',iskip
c        pause
        do k=1,iskip
          call ReadXDLine(ipar,line)
        enddo
c---  cartesian coordinates
       call rabs(0,xf(i,1),xf(i,2),xf(i,3),xc(i,1),xc(i,2),xc(i,3))
c---   goto the next atom
      enddo
c.....skip kappa sets (we don't need them)
      do i=1,nz
        call ReadXDLine(ipar,line)
      enddo
c.....read extinction parameters, out and sc
      read(ipar,*) (fvar(k),k=2,8)
      read(ipar,*)  fvar(9)
      read(ipar,*)  fvar(1)
c.....calculate the number of atom types
      ntyp=1
      iatype(1)=itype(1)      
      do 500 i=2,ia
        do j=1,ntyp
          if(itype(i).eq.iatype(j)) goto 500
        enddo
        ntyp=ntyp+1
        iatype(ntyp)=itype(i)
500   continue
      call PrAtTyp
c.....OK all done  
      call CheckOcc
      call PrintAtoms(nout)
      close(ipar)
      return
c
691   call Err( 'ReadXD', 
     *  'Error opening XD master file'//mas(:ils(mas)) )
692   call Err( 'ReadXD',
     * 'Error opening XD parameter file '//par(:ils(par)) )
693   call Err( 'ReadXD', 'Error reading XD master file')
694   call Err( 'ReadXD', 'Error reading XD parameter file')
c      
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c reads 1 line from XD file 'ixd' (skipping comments), converts it to 
c upper case and returns it to the calling subroutine
      subroutine ReadXDLine(ixd,line)
      character line*(*)
      line = ' '
119   read(ixd,'(a)',end=694,err=694) line
      if(line(1:1).eq.'!') goto 119
      call UpCase(line)
      return
694   call Err('ReadXD','Error reading XD parameter file')
      end 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c READ SHELX FILE
c
      subroutine ReadShelx(inp)
      include 'lsdb.inc'
      parameter( maxchar = 40 ) ! max number of character variables in the line
c
      character inp*60, line*255, dline*200, ans*1, alat(7)*1, line8*225
      character line61*50, line62*50
      dimension istar(maxchar)
c
      character vchar(maxchar)*8    
      common /line_extract/ nch,istar,vchar
c
      logical lrs
      logical resi11
      character*17, dimension(nat,1):: RESI1
      character*10, dimension(nat,1):: RESI2
      common /resi1/ RESI1, RESI2, resi11
      common /renameSHELXatoms/ lrs, irs
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
      common /files/ imas, nout, ixp, idb, icon     
      common /at_types/ ntyp, iatype(nat)
      data alat /'P','I','R','F','A','B','C'/
      common /shelx1/ fvar(200), wght(6), ifv
      common /ucell/ uc(6)
      dimension jsf(18)
      CHARACTER*2 SYMBAT(92)
      common /atsym/ symbat
      dimension isfn(nat), jtf(nat)
      character resi*10, atn*10
      character symlist(48)*30
c
90    format(100a)
      call izero(isfn,nat)
      call izero(jtf,nat)
      line= ' '
      dline= ' '
      line8= ' '
      do i=1,nat
        atom(i)= ' '
        asym(i)= ' '
      enddo
      ans=' '
      ntyp = 0
c
c690   write(*,'(a,$)') 
c     *' Enter SHELX input filename [leu-1.shelx]: '
c     *' Enter SHELX input filename [MCBENZ.res]: '
c     *' Enter SHELX input filename [asparm.shelx]: '
c      read(*,'(a)') inp
c      if(inp.eq.' ') inp='leu-1.shelx'
c      if(inp.eq.' ') inp='MCBENZ.res'
c      if(inp.eq.' ') inp='asparm.shelx'
c
      ish=57
      open(ish,file=inp,form='formatted',status='old',err=691)
      write(nout,'(/a,a)') ' Reading SHELX file : ',inp 
      write(imas,90) ('!',i=1,79)
      write(imas,90) '!                   CREATED BY PROGRAM << LSDB >>'
      write(imas,90) ('!',i=1,79)
      write(imas,90) 'TITLE TEST'   
      icent=0
c--- 
c---  READ "CELL" AND "SFAC" INSTRUCTIONS
c---

!! to comment
      nsym=1
      resi11 = .false.
3     read(ish,'(a)',end=4) line
        call UpCase(line)
        if(line(1:4).eq.'CELL') then
          read(line(5:100),*) rlam, uc
          write(imas,'(a,6f12.4)') 'CELL ',uc
          write(imas,'(a,f12.6)') 'WAVE ',rlam
          write(444,'(a4,3f10.4,3f9.4,1x,f8.5)') 'CELL ',uc, rlam
          write(444,'(A)') ''
!        else if(line(1:4).eq.'ZERR') then
!          read(line(10:100),*) zerr 
!          write(444,'(a,6f12.4)') 'ZERR ',zerr
        else if(line(1:4).eq.'LATT') then
c          call ReadChar(line,2,i1,i2)
          read(line(5:50),*) latt
          if(latt.ge.1) icent=1
          latt=abs(latt)
          if(icent.eq.1) write(imas,90) 'LATT C ',alat(latt)
          if(icent.eq.0) write(imas,90) 'LATT A ',alat(latt)
        else if(line(1:4).eq.'SYMM') then
          write(imas,'(a)')  line(:ils(line))
!! to comment          
          symlist(nsym)=line(5:ils(line))          
          nsym=nsym+1
        else if(line(1:4).eq.'SFAC') then
          call ReadLine(line)
c         there can be several SFAC instructions :
          do i = ntyp + 1, ntyp + nch - 1
            iatype(i) = iAtNum( vchar(i+1-ntyp) )
          enddo
          ntyp = ntyp + nch - 1
        else if(line(1:4).eq.'FVAR') then
          call DNumLine(line,inum)
          inum = inum - 1
          read(line(5:),*) (fvar(i), i=ifv+1,ifv+inum)
          ifv = ifv+inum
        else if(line(1:4).eq.'WGHT') then
          call DNumLine(line,inum)
          read(line(5:),*) (wght(i), i=1,inum-1)
        else if(line(1:4).eq.'RESI') then
          resi11=.true.
          k=0
221       read(line(5:50),*) nres, resi
         do i=1, nat  
         read(ish,'(a)',end=4) line8
        if(line8(1:4).eq.'RESI') then
        line=line8
        goto 221
        endif
        if(line8(1:4).eq.'HKLF') goto 222
        if(line8(1:4).eq.'END ') goto 222
        if(line8(1:4).eq.'REM ') goto 222
        if(line8(1:4).eq.'PART') goto 222
        if(line8(1:4).eq.'AFIX') goto 222
        if((line8(1:2).ne.'  ').and.(line8(1:4).ne.'! ').and.
     *  (line8(1:3).ne.'REM').and.(line8(1:4).ne.'PART').and.
     *  (line8(1:4).ne.'AFIX')) then
        k=k+1
        read(line8(1:10),*) atn
        write(line61,'(a6,1x,i3,2x,a5)') atn, nres, resi
        write(line62,'(a6,1x,i3)')atn, nres
        RESI1(k,1)=line61
        RESI2(k,1)=line62
        endif
         enddo 
222     continue      
        endif
      goto 3         
c---
4     continue
c      
c      write(*,*) ' ifv = ',ifv
c      stop
c
c     write something to xd_lsdb.mas
c
      write(imas,90) ('!',i=1,79)
      write(imas,90) '!'
      write(imas,90) '    MODULE *XDLSM'
      write(imas,90) '!'
      write(imas,90) 'SELECT  *model 4 0 0 0 based_on F test verbose 3'
      write(imas,90) 
     * 'SELECT  cycle 15 dampk 0.6 cmin 0.6 cmax 1. eigcut 1.d-09'
      write(imas,90) 'SAVE  deriv lsqmat cormat'
      write(imas,90) 'SOLVE  *inv diag *cond'
      write(imas,90) '!',('-',i=1,78)
      write(imas,90) 'SCAT CORE SPHV DEFV   1S  2S  3S  4S  2P  3P',
     * '  4P  3D  4D  4F  5S  5P  6S  6P  5D  7S  6D  5F  DELF''',
     * '   DELF''''  NSCTL'
      do i=1,ntyp
        call GetSFT(iatype(i),jsf)
        write(imas,41) symbat(iatype(i)),jsf, (0.00D00,j=1,3)
      enddo
41    format(a2,3x,'CHFW CHFW CSZD',1x,18i4,2f8.4,f7.3)
      write(imas,90) 'END SCAT'
      write(imas,90) '!',('-',i=1,78)
c
      write(nout,'(/a,6f10.3)') ' Cell parameters: ',uc
      write(nout,'(/2a)') ' Lattice type : ',alat(latt)
      if(icent.eq.1) write(nout,42) 'CentroSymmetric'
      if(icent.eq.0) write(nout,42) 'Acentric'
42    format(' Lattice is ',a)
c---  symmetry operations
      write(nout,'(/a)') ' Symmetry operations : '
      call GetSym(ish,nout)
      if(abs(latt).ne.1) call NonP(nout,alat(latt))
      if(icent.eq.1) call AddCentSymm(nout)
      write(nout,'(/a)') 
     * ' Final set of symmetry operations in matrix form:'     
      call PrSymOp(nout)
c---
      call PrAtTyp
c---
      call TranMat
      write(nout,345) 
      write(   *,345) 
345   format(/4x,'NO.  ATOM     IZ IU SF    XF       YF       ZF',
     *           '    OCCUP  U''s -->')      

!! to comment
      if(icent.eq.1) then
        write(444,'(A,I3,1X,A,A)') 'SYMM ',nsym,alat(latt),' CENTRO'
      else
        write(444,'(A,I3,1X,A)') 'SYMM ',nsym,alat(latt)        
      end if
      write(444,'(A)') '  X, Y, Z'
      do isym=1,nsym-1
        write(444,'(A)') symlist(isym)
      end do
      write(444,'(A)')
      write(444,'(A)') 'SCALE  1  1.0'
      write(444,'(A)')
      write(444,'(A)') 'FMULT  1.0'
      write(444,'(A)')
      write(444,'(A)') 'UOVER  0.0'
      write(444,'(A)')
      write(444,'(A)') 'SOLVT  0.0  50.0'
      write(444,'(A)')
      write(444,'(A)') 'EXTIN  0.0  GAUSSIAN  ISOT  TYP1'
      write(444,'(A)')

c---
c---  READING ATOM DATA FROM SHELX FILE
c---
      rewind(ish)
      ia=0 	
5     call RdShLine(ish,line,iend,iflag)
       if(iend.eq.1) goto 6
       dline = ' '
       dline(1:100)=line(1:100)
       if(iflag.eq.1) then
         call RdShLine(ish,line,iend,iflag)
         dline(101:200)=line(1:100)
       endif
c       write(*,'(3a)') '[',dline,']'
       ia = ia + 1   
c       write(*,*) ' ia = ',ia
       if(ia.gt.nat) call Err('ReadShelx','Too many atoms')
c---   atomic number
       call ReadChar(dline,2,i1,i2)
       read(dline(i1-1:),*) idum
       itype(ia)=iatype(idum)
c---   count number of each of atom types
       isfn(idum) = isfn(idum) + 1
       ldum = isfn(idum)       
c       write(*,*) 'idum = ',idum
c---   default multipole expansion lmax=4
       lmx(ia)=4
c---   for H-atoms - lmax=2       
       if(itype(ia).eq.1) lmx(ia)=2     
       isft(ia) = idum
c---   decide what to do with atom name       
       call ShAtNam(dline,atom(ia),iatype(idum),ia,ldum)
c---   lets figure out how many variables in dline
       call DNumLine(dline,inum) 
c---   NEW: process each variable in 'dline' in order to take of
c---   FVAR assignments
       call ProcessFVAR(dline,inum,ia)       
c---   isotropic   (there is no atom name in dline anymore)
       if( inum.eq.6 ) then              
          iut(ia)=1  ! isotropic TF of i-th atom
          read(dline,*,end=6) idum, 
     *     (xf(ia,j),j=1,3), occ(ia), u(ia,1)
           call ProcessU(ia,jtf)
c---   anisotropic
       else if( inum.eq.11 ) then
          iut(ia)=2  ! anisotropic TF of i-th atom
          read(dline,*,end=6) idum, 
     *    (xf(ia,j),j=1,3), occ(ia), (u(ia,j),j=1,3), (u(ia,j),j=6,4,-1)
c dbg only
          call CalcUeq(ia,Ueq)
c          write(*,400) ia,Ueq, (u(ia,j),j=1,3),(u(ia,j),j=6,4,-1)
c400       format(' i, Ueq(i) = ',i4,f12.6,3x,6f12.6)
c---   something is wrong - abort
       else 
          write(*,'(/a,a)') ' Error reading atom: ', atom(ia)(1:4) 
          write(*,'(a,a,a/)') ' Line <<',dline(1:ils(dline)),'>>' 
          stop
       endif
c---   Ok
       jt = 1
       if( iut(ia).eq.2 ) jt = 6
       write(   *,869) ia,atom(ia), itype(ia), iut(ia), isft(ia),
     *                 (xf(ia,j),j=1,3), occ(ia), (u(ia,k),k=1,jt)
       write(nout,869) ia,atom(ia), itype(ia), iut(ia), isft(ia),
     *                 (xf(ia,j),j=1,3), occ(ia), (u(ia,k),k=1,jt)
869    format( i6, 3x, a, 3i3, 3f9.5, f6.3, 6f7.3 )
c---   cartesian coordinates
       call rabs( 0, xf(ia,1), xf(ia,2), xf(ia,3),
     *               xc(ia,1), xc(ia,2), xc(ia,3) )
      goto 5
c---
c---  OK, ALL ATOMS READ FROM SHELX FILE
c---
6     continue 
      write(*,'(/a,i5)') ' Total atoms read from SHELX file = ',ia	
      call CheckOcc
      call PrintAtoms(nout)
      close(ish)
      return
691   call Err('ReadShelx','Error opening SHELX file...')
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c very crude: scat factor table for atom with atomic number 'iat' 
c when working w/shelx file
      subroutine GetSFT(iat,jsf)
      dimension jsf(18)
      call izero(jsf,18)
c     1s (1)
      if( iat.eq.  1) jsf( 1)=-1
      if( iat.eq.  2) jsf( 1)=-2
      if( iat.gt.  2) jsf( 1)= 2
c     2s (2) & 2p (5)
      if( iat.eq.  3) jsf( 2)= -1
      if( iat.ge.  4) jsf( 2)= -2
      if( iat.ge.5. and. iat.le.10 ) jsf(5) = -( iat-4 )
      if( iat.gt. 10) jsf( 2)=  2
      if( iat.gt. 10) jsf( 5)=  6
c     3s (3) & 3p (6)
      if( iat.eq. 11) jsf( 3)= -1
      if( iat.ge. 12) jsf( 3)= -2
      if( iat.ge.13. and. iat.le.18 ) jsf(6) = -( iat-12 )
      if( iat.gt. 18) jsf( 3)=  2
      if( iat.gt. 18) jsf( 6)=  6
c     4s (4), 3d (8), 4p (7)
      if( iat.eq. 19) jsf( 4)= -1
      if( iat.ge. 20) jsf( 4)= -2
      if( iat.ge.21. and. iat.le.30 ) jsf(8) = -( iat-20 )
      if( iat.ge.31 ) jsf(8) = -10
      if( iat.ge.31. and. iat.le.36 ) jsf(7) = -( iat-30 )
      if( iat.gt. 36) jsf( 4)=  2
      if( iat.gt. 36) jsf( 8)= 10
      if( iat.gt. 36) jsf( 7)=  6
c     5s (11), 4d(9), 5p(12)
      if( iat.eq. 53) jsf( 9)= 10
      if( iat.eq. 53) jsf( 11)= -2
      if( iat.eq. 53) jsf( 12)= -5
c
      if(iat.gt.36 .and. iat.ne.53) 
     * call Err('GetSFT','atomic number not implemented yet')
c               
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c reads line from SHELX file, checks first 4 characters against all 
c SHELX keywords, and if not a single keyword was found then it 
c returns the line back to main subroutine
c it also checks if the last symbol in line is '=' character
c iflag=0 - no '=', iflag = 1 if there is a '='
c if iflag=1 then '=' is replaced by ' '
c iend = 1 if reached the end of the file (otherwise iend=0)
      subroutine RdShLine(ii,line,iend,iflag)
      include 'lsdb.inc'
      character line*(*)
c     shelx keywords
      parameter( nkey=73 )
      character*4 skey(nkey)
      data skey/
     *'TITL','CELL','ZERR','LATT','SYMM','SFAC','DISP','UNIT','LAUE',
     *'REM ','MORE','TIME','FLAT','HKLF','OMIT','SHEL','BASF','TWIN',
     *'EXTI','SWAT','HOPE','MERG','SPEC','RESI','MOVE','ANIS','AFIX',
     *'HFIX','FRAG','FEND','EXYZ','EADP','EQIV','OMIT','CONN','PART',
     *'BIND','FREE','DFIX','DANG','BUMP','SAME','SADI','CHIV','FALT',
     *'DELU','SIMU','DEFS','ISOR','NCSY','SUMP','L.S.','CGLS','BLOC',
     *'DAMP','STIR','WGHT','FVAR','BOND','CONF','MPLA','RTAB','HTAB',
     *'LIST','ACTA','SIZE','TEMP','WPDB','FMAP','GRID','PLAN','MOLE',
     *'END '/   
c
      iend  = 0
      iflag = 0
c      
400   line = ' '
      read(ii,'(a)',end=700) line
c     skip blank lines and "our" comments     
      if(line.eq.'    '.or.line(1:1).eq.'!') goto 400
      call UpCase(line)
c     end of file 
      if(line(1:4).eq.'END ') goto 700
c         
      js=ils(line)
      do i=1,nkey
        if(skey(i)(1:4).eq.line(1:4)) then
          if(line(js:js).eq.'=') read(ii,'(a)',end=700) line
          goto 400
        endif  
      enddo        
      if(line(js:js).eq.'=')  then
        iflag=1
        line(js:js)=' '
      endif
      return
c.....reached end of the file
700   continue
      iend=1
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c determines atom name from SHELX record line
c ia - atomic number
      subroutine ShAtNam(line,anam,iat,isn,ldum)
      character line*(*), anam*(*), ajunk*10
      CHARACTER*2 SYMBAT(92)
      common /atsym/ symbat
      logical lrs
      common /renameSHELXatoms/ lrs, irs
c
      anam = ' '
      call ReadChar(line,1,i1,i2)
      if(.not.lrs) goto 20
c
c     change atom name :
c
c.....straightforward atom counting, i.e. C1 H2 O3 H4 ...
      if( irs.eq.0 ) then
        write(ajunk,10) symbat(iat),isn
c.....count atoms in each atom type, i.e. C1 H1 O1 H2 ...
      else    
        write(ajunk,10) symbat(iat),ldum
      endif
c
10    format(a2,'(',i6,')')
      call RmBlanks(ajunk)
      anam = ajunk(:ils(ajunk))
      line(1:i2)=' '
      return
c
c     keep original SHELX name, but add "(" and ")"
c
20    continue
      ij=2
      if(symbat(iat)(2:2).eq.' ') ij=1
      write(anam,'(4a)') line(1:ij),'(',line(ij+1:i2),')'
      line(1:i2)=' '
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c takes care of atomic parameters being defined using FVAR
      subroutine ProcessFVAR(dline,inum,ia)
      character dline*(*)
      dimension dum(100)
      common /shelx1/ fvar(200), wght(6), ifv
      data zero, ten / 0.00D00, 10.00D00 /
c     read variables to array dummy dum
      call rzero(dum,100)
      read(dline,*) (dum(i),i=1,inum)
c     process each variable      
      do i=2,inum
c       write(*,*) ' i, dum(i) = ',i,dum(i)
c      leave as is ( -10 < par < 10.0 )
       if( dum(i).gt.-ten .and. dum(i).lt.ten ) then 
c      just fixed ( -20 < par <= -10.0 .and. 10 <= par < 20.0 )
       else if( (dum(i).gt.-20.000.and. dum(i).le.-ten ) .or. 
     *         ( dum(i).ge.ten    .and. dum(i).lt.20.00D00 ) ) then
         dum(i) = dum(i) - ten
c      FVAR
       else
         m = nint(dum(i)/ten)
         iznak = 1
         if( dum(i).lt.zero ) iznak = -1
         m = abs(m)
         p = mod( dum(i), ten )
         if( iznak.eq. 1 ) then
           dum(i) = p*fvar(m)
         endif
         if( iznak.eq.-1 ) then
           dum(i) = p*(fvar(m)-1.00D00)
         endif
       endif 
      enddo
c     rewrite 'dline' with new parameters
      dline = ' '
      write(dline,'(5x,i4,100f12.6)') nint(dum(1)), (dum(i),i=2,inum)     
      return
      end      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c takes care of Uiso of k-th atom being defined relative to the Ueq 
c of last normal atom
      subroutine ProcessU(k,jtf)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      dimension jtf(nat)
c
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
c
      uia = u(k,1)
c     Uiso fixed at Ueq of the last normal atom (i.e. with jtf=0)
      if( uia.ge.-5.00D00 .and. uia.le.-0.5D00 ) then
        jtf(k) = 1
        do i=k-1,1,-1
          if( jtf(i).eq.0 ) then
            call CalcUeq(i,Ueq)
            u(k,1) = abs(uia)*Ueq
            return
          endif
        enddo
        stop ' ProcessU :: Could not find parent atom '
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for j-th atom calculate Ueq 
c R.X. Fischer & E. Tillmanns, Acta Cryst. (1988). C44, 775-776
      subroutine CalcUeq(j,Ueq)
      include 'lsdb.inc'
      dimension ud(3,3), rc(6)
      double precision ev(3), work(3), u1(3,3)
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
      common /ucell/ uc(6)
c
      if( iut(j).ne.2 ) then
        Ueq = u(j,1) 
        return
      endif  
c      
      ud(1,1) = u(j,1)
      ud(2,2) = u(j,2)
      ud(3,3) = u(j,3)
      ud(1,2) = u(j,4)
      ud(1,3) = u(j,5)
      ud(2,3) = u(j,6)
      ud(2,1) = ud(1,2)
      ud(3,1) = ud(1,3)
      ud(3,2) = ud(2,3)
      call r_cell(uc,rc)
      Ueq = 0.00D00
      do i =1,3
        Ueq = Ueq + ud(i,i) * ( uc(i)*rc(i) )**2
      enddo
      term1 = 2.00D00*ud(1,2)*uc(1)*uc(2)*rc(1)*rc(2)*cosd(uc(6))
      term2 = 2.00D00*ud(1,3)*uc(1)*uc(3)*rc(1)*rc(3)*cosd(uc(5))
      term3 = 2.00D00*ud(2,3)*uc(2)*uc(3)*rc(2)*rc(3)*cosd(uc(4))
      Ueq = ( Ueq + term1 + term2 + term3 ) / 3.00D00
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c converts real cell parameters (cell) to reciprocal (rcell)
c Angstroms & Degrees
	subroutine r_cell(cell,rcell)
	dimension cell(6)      ! real cell parameters
	dimension rcell(6)     ! reciprocal cell parameters
c
	vol = cell(1)*cell(2)*cell(3)*sqrt( 1 - cosd(cell(4))**2 - 
     *      cosd(cell(5))**2 - cosd(cell(6))**2 + 
     *      2*cosd(cell(4))*cosd(cell(5))*cosd(cell(6)) )
	rcell(1) = cell(2)*cell(3)*sind(cell(4))/vol
	rcell(2) = cell(1)*cell(3)*sind(cell(5))/vol
	rcell(3) = cell(1)*cell(2)*sind(cell(6))/vol
	dum1 = vol / (cell(1)*cell(2)*cell(3))
      rcell(4)  = asind( dum1/( sind(cell(5))*sind(cell(6)) ) )
      rcell(5) = asind( dum1/( sind(cell(4))*sind(cell(6)) ) )
      rcell(6) = asind( dum1/( sind(cell(4))*sind(cell(5)) ) )	
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c READ PLATON ANGSTROM FILE
c
      subroutine ReadPlaton(out)
      include 'lsdb.inc'
c      parameter( a01=0.529177249d00 )
      parameter( maxchar = 40 ) ! max number of character variables in the line
c
      character out*60, line*255
      character atom(nat)*8, asym(nat)*5
      character vchar(maxchar)*8
      integer istar(maxchar)
      common /line_extract/ nch,istar,vchar
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_types/ ntyp, iatype(nat)
      common /files/ imas, nout, ixp, idb, icon
      data zero /0.00D00/
c
c      out=' '
      line=' '
      itop=1
c691   write(*,'(a,$)') 
c     *' Enter name of PLATON file in ANGSTROM format [junk.ang] : '
c      read(*,'(a)') out
c      if(out.eq.' ') out='junk.ang'
      open(itop,file=out,form='formatted',status='old',err=691)
      write(nout,'(//a,a)') ' Reading PLATON file : ',out     
c---
c---  read file
c---
      i = 0
2     read(itop,99,end=6,err=6) line
        call ReadLine(line)
        if(vchar(1).ne.'ATOM') goto 2
        i = i +1  
        if(i.gt.nat) call Err('ReadPlaton','Too many atoms')
        occ(i) = 1.00D00
        atom(i)=vchar(2)
        read(line,*) (xc(i,k),k=1,3)           
        iatn=iAtNum(atom(i)(1:8))    
        itype(i)=iatn
c---    default multipole expansion lmax=4
        lmx(i)=4
c---    for H-atoms - lmax=2       
        if(itype(i).eq.1) lmx(i)=2     
c---
        if(i.eq.1) then
          ntyp=1
          iatype(ntyp)=iatn
        else
          iflag=0
          do j=1,ntyp
            if(iatn.eq.iatype(j)) iflag=1  
          enddo
          if(iflag.eq.0) then
            ntyp=ntyp+1
            iatype(ntyp)=iatn
          endif
        endif            
c          read(line(25:51),'(3f9.4)') (xc(i,k),k=1,3) 
c          write(*,99) 'line(25:51)=[',line(25:51),']'      
        do k=1,3
          xf(i,k)=zero
        enddo
        goto 2
c---
c---  OK, ALL ATOMS READ FROM 4TOPOND FILE
c---
6     continue 
      ia = i
      close(itop)
      call PrAtTyp
c      write(*,'(/a,i4)') ' Number of atom types = ',ntyp
c      do i=1,ntyp
c        write(*,429) i, iatype(i)   
c      enddo
c429   format(3x,'atom type = ',i2,3x,' atomic number = ',i3)     
c
      write(   *,328) ia 
      write(nout,328) ia 
328   format(/' Total atoms read from PLATON file = ',i8)
c
      call PrintAtoms(nout)
99    format(10a)
      return
691   call Err('ReadPlaton','Error opening PLATON file')
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c  READS AIM2TAB OUTPUT FILE
c
      subroutine ReadAim2Tab(out)
      include 'lsdb.inc'
      parameter( a01=0.529177249d00 )
c
      character out*60, line*255
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_types/ ntyp, iatype(nat)
      common /files/ imas, nout, ixp, idb, icon
      data zero /0.00D00/
c
c      out=' '
      line=' '
      itop=1
c691   write(*,'(a,$)') 
c     *' Enter AIM2TAB output filename [aim2tab.out] : '
c      read(*,'(a)') out
c      if(out.eq.' ') out='aim2tab.out'
      open(itop,file=out,form='formatted',status='old',err=691)
c---
c---  read file
c---
2     read(itop,99,end=6,err=6) line
        if(line(2:23).eq.'NUMBER OF ATOMS READ =') then
          read(line(24:30),*) ia 
          if(ia.gt.nat) call Err('ReadAim2Tab','Too many atoms')
        else if(line(2:23).eq.'NO. ATOM       AT.MASS') then
          read(itop,99) line
          do i=1,ia
           read(itop,'(5x,a8,11x,3f9.4)') atom(i)(1:8), (xc(i,k),k=1,3)           
           occ(i) = 1.00D00
c           write(*,*) (xc(i,k),k=1,3)
c           atom(i)(1:8)=line(6:13)
           iatn = iAtNum(atom(i)(1:8))    
           itype(i)=iatn
c---       default multipole expansion lmax=4
           lmx(i)=4
c---       for H-atoms - lmax=2       
           if(itype(i).eq.1) lmx(i)=2     
c---
           if(i.eq.1) then
             ntyp=1
             iatype(ntyp)=iatn
           else
             iflag=0
             do j=1,ntyp
               if(iatn.eq.iatype(j)) iflag=1  
             enddo
             if(iflag.eq.0) then
               ntyp=ntyp+1
               iatype(ntyp)=iatn
             endif
           endif            
c           read(line(25:51),'(3f9.4)') (xc(i,k),k=1,3) 
c           write(*,99) 'line(25:51)=[',line(25:51),']'
           do k=1,3
             xf(i,k)=zero
             xc(i,k)=xc(i,k)*a01
           enddo
          enddo 
        endif
      goto 2
c---
c---  OK, ALL ATOMS READ FROM 4TOPOND FILE
c---
6     continue 
      close(itop)
      write(*,'(/a,i4)') ' Number of atom types = ',ntyp
      do i=1,ntyp
        write(*,429) i, iatype(i)   
      enddo
429   format(3x,'atom type = ',i2,3x,' atomic number = ',i3)     
      write(*,'(/a,i5)') 
     * ' Total atoms read from AIM2TAB output file = ',ia	
      call PrintAtoms(nout)
99    format(10a)
      return
691   call Err('ReadAim2Tab','Error opening AIM2TAB output file...')
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c reading and decoding SYMM cards from input file (taken from XD)
c
c***********************************************************************
c
C---- Reads the SYMM entries in the Master File (unit fileno),
C---- stores the original set of ns given symmetry operations in arrays
C---- fs (rotation) and ts (translation),
      Subroutine GetSym(fileno, nout)
      Logical eof, lrew
      Integer   fileno
      Dimension t(3,4)
      Character*255 rec
      common /symm/ ns, fs(3,3,48), ts(3,48)
      data zero /0.00D00/
C
C---- Default is centrosymmetric, primitive lattice
C
      Do i=1,3
         t(i,1)=zero
      End Do
C
C---- Decode lines beginning with SYM(M) in the Master File
C---- SYMM X,Y,Z is not explicitly required and will be ignored if given.
C---- It is always operation number 1.
C
      Do i=1,3
         ts(i,1)   = zero
         fs(i,i,1) = 1.
         Do j=1,3
            If (i .Ne. j) fs(i,j,1) = zero
         Enddo
      Enddo
      ns = 1
      Rewind fileno
      lrew = .False.
 100  Continue
      Call SRCHL(fileno, 'SYM', lrew, rec, eof)
      If (.Not. eof) Then
         Call SYMCRD(rec, ts, fs, ns, ierr)
         If (ierr .Ne. 0) Then
            Write(nout,9001) ns
 9001       Format(' * * * Symmetry card no.', I3, ' illegal * * *')
            ns=ns-1
         Endif
         Goto 100
      Endif
      Do i=1,ns
         Call SYMOUT(nout, i, fs, ts)
      Enddo
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine SYMOUT(lfc, is, rs, ts)
      Dimension rs(3,3,*), ts(3,*)

      Character*80 line
      Integer start

      start=1
      line = ' '
      Do 100 J=1,3
         Call OPOUT(line, start, rs(1,j,is), rs(2,j,is), rs(3,j,is),
     U        ts(j,is))
         If (j .Eq.3) Goto 100
         line(start:start)=','
         start=start+1
 100  Continue
      start=start-1
      Write(lfc,9010) is, line(:start)
 9010 Format(3x,'SYMM', I3, 5X, A)
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine OPOUT(line, i, a1, a2, a3, t)
      Character*80 line

      Dimension out(3)
      Character*1 tag(3)
      Data tag/'X','Y','Z'/

      out(1) = a1
      out(2) = a2
      out(3) = a3
      If (t .Ne. 0.) Then
         k=Int(12*t+0.5)
         Do 100 l=1,6
            If (l .Ne. 3) Then
               If (l .Eq. 2) Then
                  j = 4
               Else
                  j = 7-l
               Endif
               If (Mod(k,j) .Eq. 0) Then
                  k =  k/j
                  j = 12/j
                  Write(line(i:i+3),1010) k, j
 1010             Format( I2, '/', I1)
                  i=i+4
                  Goto 190
               Endif
            Endif
 100     Continue
      Endif
 190  Continue

      Do 200 J=1,3
         If (out(j) .Eq. 0.) Goto 200
         If (out(j) .Gt. 0.) Then
            line(i:i+2) = ' + '
         Else
            line(i:i+2) = ' - '
         Endif
         i = i+3
         a = Abs(out(j))
         If (a .Ne. 1.) Then
            Write(line(i:i+2),1020) a
 1020       Format(F3.1)
            i=i+3
         Endif
         line(i:i) = tag(j)
         i=i+1
 200  Continue
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine SYMCRD(line, ts, rs, ns, ierr)
      Character*(*) line
      Dimension rs(3,3,*), ts(3,*)

      Character*80 sym(3)
      Character*32 op(12)
      Character*4  lab
      Dimension lsym(3), lop(12)
      
      i0 = 1
      lab = line(:4)
      Call UPCASE(lab)
      If ( lab(:3) .Eq. 'SYM' ) i0 = 4
      If ( lab(:4) .Eq. 'SYMM') i0 = 5

      ms = 3
      Call SPLIT(line(i0:), ';', ms, sym, lsym, 3)
      Do 100 i=1,ms
         ns = ns+1
         mo = 12
C      Print *, 'symop: ', sym(i), lsym(i)
         If (Index(sym(i),',') .Eq. 0) Then
            Call SPLIT(sym(i)(:lsym(i)), ' ', mo, op, lop, 1)
         Else
            Call SPLIT(sym(i)(:lsym(i)), ',', mo, op, lop, 1)
         Endif
         If (mo .Eq. 12) Then
            l = 0
            Do j=1,3
               l=l+1
               Read(op(l)(:lop(l)),*,Err=900) v
               ts(j,ns) = v
               Do k=1,3
                  l=l+1
                  Read(op(l)(:lop(l)),*,Err=900) v
                  rs(k,j,ns) = v
               Enddo
            Enddo
         Else If (mo .Ne. 3)  Then
            Goto 900
         Else
            Do j=1,mo
C               Print *, '   op: ', op(j)(:lop(j)), lop(j)
               Call SYMOP(op(j)(:lop(j)), ts(j,ns), rs(1,j,ns))
            Enddo
         Endif
C      Print 6001, ns, (ts(j,ns),j=1,3), ns, 
C     +        ((rs(k,j,ns),k=1,3),j=1,3)
C 6001 Format( 'ts[', i2, ']', 2X, 3F4.1, '   rs[', i2, ']',
C     +        3(2X, 3F4.1))
C
C     ignore identity operation
C                                    
        Do j=1,3
           If (rs(j,j,ns) .Ne. 1.) Goto 100
           If (ts(j,ns)   .Ne. 0.) Goto 100
           Do k=1,3
              If (j .Ne. k .And. rs(j,k,ns) .Ne. 0.) Goto 100
           Enddo
        Enddo
        ns=ns-1

 100  Continue

      ierr = 0
C      Print *, ns
      Return

 900  Continue
      ierr = 1021
      Return
      End
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c decode a single symmetry operator. first the rotation part is 
c determined and removed, then - translation. i think it's a much better
c way of doing it, compared to original XD subroutine
      Subroutine SYMOP(op, t, r)
      Character op*(*), top*32, c*1, num(4)*1
      Dimension r(3), y(2)
c      common /files/ imas, io, ixp, idb 
      Data num /'1', '2', '3', '4'/
cdbg
c      write(io,'(/4a)') 'symop: op=[',op(1:ils(op)),']'
      top=' '
      call RmBlanks(op)
c      write(io,'(3a)')  '       op=[',op(1:ils(op)),']'
      call RemoveRot(op,r)
      id=min(len(op),len(top))
      top(1:id)=op(1:id)
C
C 3) Analyze the translational part (first as a real number, then
C    as a quotient)
C
      t = 0.
C      Print *, '(x)    t:', t, '   r:', r
      If (top .Eq. ' ') Return
      value = 0.
      If (Index(top, '/') .Ne. 0) Goto 100
      Read(top,*,err=100,end=100) value
      t = value
C      Print *, '(f)    t:', t, '   r:', r
      Return

 100  Continue
      indx=1
      y(1)=0.
      y(2)=1.
      sgn = 1.
      Do 110 l=1,LENSTR(top)
         c = top(l:l)
         If (c .Eq. ' ') Goto 110
         If (c .Eq. '/') indx=2
         Do k=1,4
            If (c .Eq. num(k)) y(indx)=sgn*Float(k)
         Enddo
         t=y(1)/y(2)
         sgn=1.
         If (c .eq. '-') sgn=-1.
  110  Continue

cdbg
c      write(io,*) '(q)    t:', t, '   r:', r
cenddbg
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determines the rotation part from symmetry element "op", puts it in
c numerical array "r(3)" and "removes" the rotation part from "op".
c thus, only a translation part is left in "op"
      subroutine RemoveRot(op,r)
      character axp(3)*2, axn(3)*2, ax(3)*1, op*(*)
      Dimension r(3)
      data axp /'+X','+Y','+Z'/, axn /'-X','-Y','-Z'/, ax  /'X','Y','Z'/
      call rzero(r,3)
c.....loop over possibilities +?, -?, ?
      do j=1,3
c.......loop over axes X, Y, Z      
        do i=1,3
          if(j.eq.1) id=index(op,axp(i))
          if(j.eq.2) id=index(op,axn(i))
          if(j.eq.3) id=index(op, ax(i))
          il=2
          if(j.eq.3) il=1  
          if(id.ne.0) then
            r(i)=1.00D00
            if(j.eq.2) r(i)=-1.00D00
            op(id:id+il-1)=' '
            call RmBlanks(op)
          endif  
        enddo
      enddo  
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine SRCHL(file, key, rewfil, line, eof)
      implicit double precision (a-h,o-z)
C
C Search for line starting with 'key' and read.
C
      Parameter (MAXKEY=32)
      Character key*(*), line*(*)
      Character*(MAXKEY) uk, ul
      Integer   file
      Logical   rewfil, eof
C
      eof = .False.
      lk = Len(key)
C      lk  = Max(4,lk)
      uk  = key
      Call UPCASE(uk(:lk))
      If (rewfil) Rewind file
      rewfil = .False.
C
 100  Continue
      Read(file,9010,End=190) line
      ul = line
      Call UPCASE(ul)
      If (uk(:lk) .eq. ul(:lk)) Return
      Goto 100
C
 190  Continue
C      Print 9020, file, key
      eof = .true.
      Return 
C
 9010 Format(A)
 9020 Format(' EOF on file ', I2, ' while searching for ', A,
     + ' (may be expected)')
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine SPLIT(string, sep, maxtok, tok, length, iflags)
c
      Integer STRIP, UPPER, LOWER
      Parameter (STRIP=1)
      Parameter (UPPER=2)
      Parameter (LOWER=4)

      Character*(*) string, tok(maxtok)
      Character*1   sep, nul, tab
      Integer maxtok, length(maxtok), flags

      Logical lstr, llow, lupp, lcmp

      maxlen = Len(tok(1))
      flags = iflags
      lstr = .False.
      llow = .False.
      lupp = .False.
      lcmp = .False.
      nul = Char(0)
      tab = Char(9)
      If (Mod(flags,STRIP+1) .Eq. STRIP) Then
         flags = flags-STRIP
         lstr = .True.
      Endif
      If (Mod(flags,UPPER+1) .Eq. UPPER) Then
         flags = flags-UPPER
         lupp = .True.
      Endif
      If (Mod(flags,LOWER+1) .Eq. LOWER) Then
         flags = flags-LOWER
         llow = .True.
      Endif
      If (lstr .And.
     +     (sep .Eq. ' ' .Or. sep .Eq. nul .Or. sep .Eq. tab)) Then
         lcmp = .True.
      Endif

      i=0
      j=1
      n=0
      l=Len(string)
 100  Continue
         i = i+1
         If (string(i:i) .Eq. sep) Then
            n=n+1
            if (i .Eq. j) Then
               If (lcmp) Then
                  n = n - 1
                  j = i + 1
               Else
                  tok(n) = ' '
                  length(n) = i - j
                  j = i+1
                  If (n .Ge. maxtok) Goto 190
               Endif
            Else
               If (i-j .Gt. maxlen) Then
                  m = j+maxlen-1
               Else
                  m = i-1
               Endif
               tok(n) = string(j:m)
               length(n) = m+1 - j
               j = i+1
               If (n .Ge. maxtok) Goto 190
            Endif
         Else If (i .Ge. l) Then
            n=n+1
            If (i+1-j .Gt. maxlen) Then
               m = j+maxlen-1
            Else
               m = i
            Endif
            tok(n) = string(j:m)
            length(n) = m+1 - j
            Goto 190
         Endif
         If (i .Lt. l) Goto 100

 190  Continue
      maxtok = n

      If (.Not. (lstr .Or. llow .Or. lupp)) Return

      Do 300 n=1,maxtok
         l = length(n)
         If (lstr) Then
            j = 0
 210        Continue
            j = j + 1
            If (j .Le. l .And. (tok(n)(j:j) .Eq. ' ' .Or.
     +           tok(n)(j:j) .Eq. nul .Or.
     +           tok(n)(j:j) .Eq. tab)) Then
               Goto 210
            Endif

            If (j .Le. l) Then
               i = l+1
 220           Continue
               i = i - 1
               If (i .Gt. j .And. (tok(n)(i:i) .Eq. ' ' .Or.
     +              tok(n)(i:i) .Eq. nul .Or.
     +              tok(n)(i:i) .Eq. tab)) Then
                  Goto 220
               Endif

               i = i - 1
               j = j - 1
               length(n) = i - j + 1
               If (j .Gt. 0) Then
                  Do 230 k=1,length(n)
                     tok(n)(k:k) = tok(n)(k+j:k+j)
 230              Continue
               Endif
               tok(n)(length(n)+1:) = ' '
            Else
               tok(n)    = ' '
               length(n) = 0
            Endif
         Endif

         k = 1
         If (lupp) Then
            Call UPCASE(tok(n)(k:))
            k = 2
         Endif
         If (llow) Then
            Call LOCASE(tok(n)(k:))
         Endif
 300  Continue

      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c adds center of symmetry to symmetry elements          
      subroutine AddCentSymm(nout)
      common /symm/ ns, fs(3,3,48), ts(3,48)
      do 551 k=1,ns
        do 551 i=1,3 
          ts(i,k+ns)=ts(i,k)
          do 551 j=1,3 
            fs(i,j,k+ns) = -fs(i,j,k)
551   continue     
      ns = ns*2
      write(nout,'(/a,i4)') 
     * ' Symmetry operations with Center of Symmetry applied: ',ns
      do k=1,ns
        write(nout,'(/a,i4)') ' SYMM ',k
        do i=1,3
          write(nout,'(3f12.6,5x,f12.6)') (fs(i,j,k),j=1,3), ts(i,k)
        enddo
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c adds "latmod" lattice centering to symmetry elements
      subroutine NonP(nout,latmod)
      character latmod*1
      dimension tc(3,4)
      common /symm/ ns, fs(3,3,48), ts(3,48)
      write(nout,'(/a)') ' Non-Primitive lattice...'
      write(nout,'(3a)') ' Lattice centring type is [ ',latmod,' ]'
      data rhalf /0.5D00/
c
      call rzero(tc,12)
c
      if(latmod.eq.'C') then
         tc(1,1) = rhalf
         tc(2,1) = rhalf 
         latt = 2
      else if(latmod.eq.'A') then
         tc(2,1) = rhalf 
         tc(3,1) = rhalf
         latt = 2
      else if(latmod.eq.'B') then
         tc(1,1) = rhalf 
         tc(3,1) = rhalf 
         latt = 2
      else if(latmod.eq.'I') then
         tc(1,1) = rhalf
         tc(2,1) = rhalf
         tc(3,1) = rhalf 
         latt = 2
      else if(latmod.eq.'F') then
         tc(1,1) = rhalf 
         tc(2,1) = rhalf
         tc(2,2) = rhalf
         tc(3,2) = rhalf
         tc(1,3) = rhalf
         tc(3,3) = rhalf
         latt = 4
      else
         write(*,*)' Lattice type ',latmod,' is not supported'
         stop
      endif
      write(nout,'(a,i1)') ' Number of lattice points per cell = ',latt
c
      j1=0
c
      do 601 i=1,latt-1
       do 601 j=1,ns
        j1=i*ns+j
cccc        write(nout,'(3(a,i4,5x))') 'i=',i,'j=',j,'j1=',j1
        do 601 k=1,3
          ts(k,j1)=ts(k,j)+tc(k,i)
          do 601 l=1,3
            fs(k,l,j1)=fs(k,l,j)
601   continue
c
      ns=j1
c---  make sure that translations are within -1...1 
      call CheckTs(ns,ts)
c
      write(nout,'(/a)')
     * ' Symmetry operations after applying the lattice centering:'
      call PrSymOp(nout)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c make sure that all given translations are within -1...1 range
      subroutine CheckTs(nj,tj)
      dimension tj(3,48)
      do 300 i=1,nj
      do 300 j=1,3
        if(tj(j,i).ge.1.d0.or.tj(j,i).le.-1.d0) 
     *                       tj(j,i) = tj(j,i) - int( tj(j,i) )         
300   continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine PrSymOp(nout)
      common /symm/ ns, fs(3,3,48), ts(3,48)
      common /verbosity/ iverbose
      write(nout,'(/a,i4)')' Number of symmetry operations = ',ns
      if( iverbose.le.1 ) return
      do k=1,ns
        write(nout,'(/a,i2)') ' SYMM ',k
        do i=1,3
          write(nout,'(3f12.6,5x,f12.6)') (fs(i,j,k),j=1,3), ts(i,k)
        enddo
      enddo
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c create matrices for CRYSTAL <--> CARTESIAN transformations 
c unit cell parameters are from common /uc/ (Angstroms and Degrees)
      subroutine TranMat
      Parameter (pi=3.141592653589793238d0,torad=pi/180.0d0)
      dimension uj(3,3)
      common /files/ imas, nout, ixp, idb, icon
      common /transf/ rm(3,3), rm1(3,3)
      common /ucell/ uc(6)
      data zero /0.00D00/
      a=uc(1)
      b=uc(2)
      c=uc(3)
      alpha=uc(4)
      beta=uc(5)
      gamma=uc(6)
      alpha=torad*alpha
      beta=torad*beta
      gamma=torad*gamma
      uj(1,1)=a
      uj(2,1)=zero
      uj(3,1)=zero
      uj(1,2)=b*cos(gamma)
      uj(2,2)=b*sin(gamma)
      uj(3,2)=zero
      uj(1,3)=c*cos(beta)
      uj(2,3)=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
      uj(3,3)=sqrt(c**2-uj(1,3)**2-uj(2,3)**2)
c      V=a*b*c*sqrt(1.d0-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+
c     U  2.d0*cos(alpha)*cos(beta)*cos(gamma))
      do i=1,3
       do j=1,3
        rm(i,j)=uj(i,j)
       enddo 
      enddo  
c
      write(nout,'(/a)') 
     *' Transformation matrix CRYSTAL --> GLOBAL CARTESIAN :'
      do i=1,3
        write(nout,'(3f15.8)') (rm(i,j),j=1,3)
      enddo
      call Inverse(rm,rm1)
      write(nout,'(/a)') 
     *' Transformation matrix GLOBAL CARTESIAN --> CRYSTAL :'
      do i=1,3
        write(nout,'(3f15.8)') (rm1(i,j),j=1,3)
      enddo
      write(nout,'(a)') '  '
      write(   *,'(a)') '  '
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
c    CALCULATES INVERSE OF 3x3 MATRIX: C=Inverse[A]  (NOT VERY SMART!!!)
      subroutine Inverse(a,c)
      dimension a(1:3,1:3),c(1:3,1:3)
      data zero, one / 0.00D00, 1.00D00 /
c
      det= a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) - 
     *     a(2,1)*a(1,2)*a(3,3) + a(2,1)*a(1,3)*a(3,2) + 
     *     a(3,1)*a(1,2)*a(2,3)  -a(3,1)*a(1,3)*a(2,2)
c
      if( det.eq.zero ) call Err('Inverse','Determinant = 0')
      det = one/det
c
      c(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
      c(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
      c(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
      c(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      do i=1,3
        do j=1,3
         c(i,j)=c(i,j)*det
        enddo 
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c test each atom that does not have an occupancy of 1.00 for special
c position and if not identified - assigned as disordered one...
c iaf(i) =  0 : normal atom
c        = +1 : in special position
c        = -1 : disordered  
      subroutine CheckOcc
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /symm/ ns, fs(3,3,48), ts(3,48)
      common /files/ imas, nout, ixp, idb, icon 
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      dimension x1(3), x2(3)
      data one /1.00D00/
c     loop over atoms
      do 2 i=1,ia
        iaf(i) = 0
        if( occ(i).eq.one ) goto 2
c---    for each symmetry operation:
        do 900 is=1,ns
c---      additional translations along each of the axis     
          do 901 itx=-2,2,1
          do 901 ity=-2,2,1
          do 901 itz=-2,2,1  
c----        skip identity translation      
             if(is.eq.1.and.itx.eq.0.and.ity.eq.0.and.itz.eq.0) goto 901
c----        Coordinates of the new symmetry-related atom 
             x1(1) = xf(i,1)*fs(1,1,is) + xf(i,2)*fs(2,1,is)
     *             + xf(i,3)*fs(3,1,is) + ts(1,is) + real(itx)
             x1(2) = xf(i,1)*fs(1,2,is) + xf(i,2)*fs(2,2,is)
     *             + xf(i,3)*fs(3,2,is) + ts(2,is) + real(ity)
             x1(3) = xf(i,1)*fs(1,3,is) + xf(i,2)*fs(2,3,is)
     *             + xf(i,3)*fs(3,3,is) + ts(3,is) + real(itz)
c----        Transform to cartesian
             call rabs(0,x1(1),x1(2),x1(3),x2(1),x2(2),x2(3))               
c----        check if the generated atom is the i-th atom itself,
c----        like when i-th atom is in special position
             r1 = Dist( xc(i,1), xc(i,2), xc(i,3), x2(1), x2(2), x2(3) )
             if( r1.lt.1.00D-03 ) then
                iaf(i) = 1
                goto 2
             endif         
901       continue
900     continue
        iaf(i) = -1 
c        
2     continue      
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c print out the number of atom types
      subroutine PrAtTyp
      include 'lsdb.inc'
      character*2 symbat(92)
      common /atsym/ symbat
      common /at_types/ ntyp, iatype(nat)
      common /files/ imas, nout, ixp, idb, icon
      write(nout,'(/a,i4)')     ' Number of atom types = ',ntyp
      do i = 1,ntyp
        write(nout,429) symbat(iatype(i)), iatype(i)             
      enddo
429   format(3x,'atom type = ',a,3x,' atomic number = ',i3)     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c printout atom list to unit 'io'
      subroutine PrintAtoms(io)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /at_uij/ u(nat,6)
      write(io,329)
      ispp = 0
      idis = 0
      iord = 0
329   format(/9x,'ATOM',6x,'------- FRACTIONAL -------',4x,
     *       '------- CARTESIAN -------  FLAG')
      do i=1,ia
         write(io,330) 
     *       i, atom(i), (xf(i,k),k=1,3), (xc(i,k),k=1,3), iaf(i) 
         if( iaf(i).eq. 0 ) iord = iord + 1
         if( iaf(i).eq. 1 ) ispp = ispp + 1
         if( iaf(i).eq.-1 ) idis = idis + 1
      enddo
      write(io,333) iord, ispp, idis
333   format(/,       
     * ' FLAG =  0 : normal atoms ( ',i5,' )',/,  
     * '      =  1 : atoms in special positions ( ',i5,' ) ',/,
     * '      = -1 : disordered atoms ( ',i5,' ) ')
330   format(i7,2x,a,1x,3f9.5,2x,3f9.4,i6)      

      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c generates complete fragment(s) for CRYSTAL systems 
c info stored in : commons /atoms/ & /cryst/ 
      subroutine GenCompleteFrags
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /symm/ ns, fs(3,3,48), ts(3,48)
      common /files/ imas, nout, ixp, idb, icon 
      common /covrad/ covrad(54)
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /tols/ tol1, tol2, tol3
c     natoms : number of atoms in COMPLETE fragment(s), including all 
c              generated symmetry-equivalents
c     isymeq : 0 - symmetry-independent; 
c              otherwise the number of symm.-indep. from which it was generated
c     nsymeq : for each symm.-indep. atom the number of symm.-equiv. atoms 
c              that were generated from it    
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      dimension x1(3), x2(3)
      common /compfrag/ ncells

      natoms = ia

1     continue
c     grow fragment atom by atom 
      do 2 i=1,ia
        itypi = itype(i)
c        write(*,'(/a,i4)') ' atom = ',i
c----   for each symmetry operation:
        do 900 is=1,ns
c---      additional translations along each of the axis     
          do 901 itx=-ncells,ncells,1
          do 901 ity=-ncells,ncells,1
          do 901 itz=-ncells,ncells,1  
c----        skip identity translation      
             if(is.eq.1.and.itx.eq.0.and.ity.eq.0.and.itz.eq.0) goto 901
c----        Coordinates of the new symmetry-related atom 
             x1(1) = xf(i,1)*fs(1,1,is) + xf(i,2)*fs(2,1,is)
     *             + xf(i,3)*fs(3,1,is) + ts(1,is) + real(itx)
             x1(2) = xf(i,1)*fs(1,2,is) + xf(i,2)*fs(2,2,is)
     *             + xf(i,3)*fs(3,2,is) + ts(2,is) + real(ity)
             x1(3) = xf(i,1)*fs(1,3,is) + xf(i,2)*fs(2,3,is)
     *             + xf(i,3)*fs(3,3,is) + ts(3,is) + real(itz)
c----        Transform to cartesian
             call rabs(0,x1(1),x1(2),x1(3),x2(1),x2(2),x2(3))               
c----        check if the generated atom is not the i-th atom itself,
c----        like when i-th atom is in special position
             r1 = dist(xc(i,1),xc(i,2),xc(i,3),x2(1),x2(2),x2(3))
             if( r1.lt.1.00D-03 ) goto 901
c----        Loop over all (even symm.-gen.) atoms 
c----        to check if this atom is already in the list
             do j=1,natoms
c----          calculate distance  
               r1 = dist(xc(j,1),xc(j,2),xc(j,3),x2(1),x2(2),x2(3))
c----          check if this atom is already in the list
               if( r1.lt.1.00D-03 ) goto 901
             enddo   
c---         Ok, this symm.-equiv. is NOT yet in the list
c---         Let's check it's connectivity          
             do 903 j=1,natoms
c----          calculate distance  
               r1 = dist(xc(j,1),xc(j,2),xc(j,3),x2(1),x2(2),x2(3))
c----          Check the sum of radii   
               itypj = itype(j)
               rdum = covrad( itypi ) + covrad( itypj ) + tol1
               if( r1.gt.rdum ) goto 903
c---           ok, we have a new symmetry-equivalent atom in the fragment
c               write(*,*) ' is,itx,ity,itz,j,r1=',is,itx,ity,itz,j,r1 
               goto 800
903          continue
901       continue
900     continue
2     continue
      goto 888

800   continue
c---  ok, we have a new symmetry-equivalent atom in the fragment!
c     add it to array of atoms and start over         
      natoms = natoms + 1
c      write(*,*) ' natoms = ',natoms
      nsymeq(i) = nsymeq(i) + 1
      isymeq(natoms) = i
c     create atom symm.-equiv. name 
      atom(natoms) = atom(i) 
c      ajunk = ' '
c      ij = index(atom(i),')')      
c      if( ij.eq.0 ) then
c        write(ajunk,'(a,i3)') atom(i), nsymeq(i)
c      else
c        write(ajunk,'(a,i3,a)') atom(i)(1:ij-1), nsymeq(i), atom(i)(ij:)
c      endif     
c      call RmBlanks(ajunk)
c      atom(natoms)(1:8) = ajunk(1:8)
c      
      itype(natoms) =  itype(i)
      occ(natoms) = occ(i)
      iaf(natoms) = iaf(i)
      do j=1,3
        xf(natoms,j) = x1(j)
        xc(natoms,j) = x2(j)
      enddo
cDBG      
c99    format(100a)
c      idum = 44      
c      open(idum,file='junk.ang')
c      write(idum,99) 'ANGSTROM'
c      do i=1,natoms
c        write(idum,45) atom(i), (xc(i,k),k=1,3)
c      enddo
c45    format('ATOM ',a10,3x,3f12.4)
c      close(idum)
c      if( natoms.eq.ia+6 ) stop
cendDBG
      goto 1
            
c---  ok, we have generated complete fragments now
888   continue 
c     do not print new list if we did not add any atoms....   
      if(ia.eq.natoms) return
c     printout list of atoms in complete fragment(s)
      write(nout,'(/1x,a)') 'List of atoms in complete fragment(s) :'
      idum = 44      
      open(idum,file='complete.ang')
      write(nout,329) 
329   format(/7x,'ATOM',8x,'------ FRACTIONAL ------',4x,
     *       '----- CARTESIAN -----',5x,'ISQ NSQ')
      do i=1,natoms
         write(nout,330) i, atom(i), (xf(i,k),k=1,3), (xc(i,k),k=1,3),
     *                   isymeq(i), nsymeq(i)       
         write(idum,45) atom(i), (xc(i,k),k=1,3)
      enddo
      write(nout,'(/a)') 
     * ' File with complete fragment(s) in PLATON format : complete.ang' 
45    format('ATOM ',a10,3x,3f12.4)
      close(idum)
330   format(i6,1x,a,1x,3f9.4,1x,3f8.3,i8,i4)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c 
c       RING SEARCH BASED ON EXPLORATION OF CONNECTIVITY TABLE
c
c***********************************************************************
c create connectivity matrix (icmat) of the "new" atom array
c the "new" atom array is built by step-by-step removal
c of atoms with only 1 connection
      subroutine Connect4Rings
      include 'lsdb.inc'
c
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /files/ imas, nout, ixp, idb, icon 
c     "new" atom array 
      common /newat/ ia1, iaref(nat)
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      character adum*12
      character symbat(92)*2
      common /atsym/ symbat
      common /verbosity/ iverbose
c
c copy data from "main" array to array1
c
      ia1 = natoms
      do i=1,ia1
        iaref(i) = i
      enddo
      write(nout,'(/80a)') ('-',i=1,80)
      write(nout,345)
      write(   *,345)
345   format(/' CREATING CONNECTIVITY TABLE AND ELIMINATING OPEN ATOMS')
      ic = -1
c
c "elimination" loop
c
1     continue 
        ic = ic + 1
        write(nout,358) ic, ia1
        write(   *,358) ic, ia1
c        write(nout,*) ' elimination cycle ',ic,ia1
        call Eliminate(ic,ia2)
        idel = ia1-ia2
        write(nout,359) idel
        write(   *,359) idel
        ia1 = ia2
      if( idel.gt.0 ) goto 1
358   format(/' **** ELIMINATION CYCLE =',i5,5x,'NUMBER OF ATOMS =',i7)
359   format(8x,'NUMBER OF ATOMS REMOVED FROM THE TABLE = ',i7)
c
c     Ok, only atoms with 2 and more connections are left
c  
      if( iverbose.ge.2 ) then
        write(nout,'(/a)') 
     *  ' CONNECTIVITY TABLE AFTER THE FINAL ELIMINATION CYCLE :'
        call WriteConn
      endif  
      open(67,file='reduced.ang',form='formatted')
      open(68,file='reduced_1.ang',form='formatted')
      write(67,'(a)') 'ANGSTROM'
      write(68,'(a)') 'ANGSTROM'
      do 903 i=1,ia1
        iat = iaref(i)
c       platon file in which we keep the labels  
        write(67,788) atom(iat), (xc(iat,k),k=1,3)
c       platon file in which we change the labels
        itypi = itype(iat)
        write(adum,'(a,i6)') symbat(itypi), i
        call RmBlanks(adum)
        write(68,788) adum, (xc(iat,k),k=1,3)
c
903   continue
788   format('ATOM ',a8,3f12.4)
      close(67)
      close(68)
c   
c     Here's another trick: in 'icmat' array we rearrange
c     atoms so that first neighbors those with highest connectivity 
c
c      call SortICMAT      
c   
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      subroutine SortICMAT
      include 'lsdb.inc'
      common /newat/ ia1, iaref(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /files/ imas, nout, ixp, idb, icon
      logical la
c
      do i=1,ia1
c        write(nout,*) ' i = ',i
        write(nout,'(/i5,5x,10i3)') i, (icmat(i,j), j=1,nba(i) )    
        write(nout,40)            ( nba(icmat(i,j)), j=1,nba(i) )    
10      la = .false.
        do j=1,nba(i)-1
          iat1 = icmat(i,j) 
          iat2 = icmat(i,j+1) 
          if( nba(iat1).lt.nba(iat2) ) then
            icmat(i,j)   = iat2
            icmat(i,j+1) = iat1
            la = .true.          
          endif   
        enddo
        if( la ) goto 10
        write(nout,40) (icmat(i,j), j=1,nba(i) )    
        write(nout,40) ( nba(icmat(i,j)), j=1,nba(i) )    
      enddo 
40    format(10x,10i3)
c
      return
      end       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c eliminates atoms with less < 2 covalent connections
c if ic = 0 - printout connectivity matrix
      subroutine Eliminate(ic,ia2)
      include 'lsdb.inc'
      common /newat/ ia1, iaref(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /covrad/ covrad(54)
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /tols/ tol1, tol2, tol3
      common /files/ imas, nout, ixp, idb, icon
      common /verbosity/ iverbose
c
c     create connectivity matrix for a given array of atoms
c
c      write(nout,*) ' subroutine elim = ',ic,ia1
      do 900 i=1,ia1
        iat = iaref(i)
        xi = xc(iat,1)
        yi = xc(iat,2)
        zi = xc(iat,3)
        itypi = itype(iat)
        nba(i) = 0
c       get rid of disorder : 
        if( iaf(iat).eq.-1 ) goto 900
c
        do 901 j=1,ia1 
          jat = iaref(j)  
c         same atom
          if( i.eq.j ) goto 901      
c         check distance
          r = Dist( xi, yi, zi, xc(jat,1), xc(jat,2), xc(jat,3) )
          rcov = covrad( itypi ) + covrad( itype(jat) ) + tol1
c          write(nout,444) i,j,iat,jat,atom(iat),atom(jat),r,rcov
          if( r.gt.rcov ) goto 901
          nba(i) = nba(i) + 1
          icmat(i,nba(i)) = j             
901     continue  
c        write(nout,'(2(ai4))') ' atom ',i,' nba = ',nba(i) 
900   continue
c
c444   format(2i4,2x,2i5,2x,a,1x,a,2x,2f10.4)
c
      if( ic.eq.0 .and. iverbose.ge.3 ) call WriteConn
c      call WriteConn
c
      ia2 = 0  
c
c     keep only atoms with 2 and more connections
c
      do 902 i=1,ia1
        if( nba(i).lt.2 ) goto 902
        ia2 = ia2 + 1
        iaref(ia2) = iaref(i)
902   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c writes out the connectivity matrix with 'iaref' pointers to atoms in
c the main list /atoms/
      subroutine WriteConn
      include 'lsdb.inc'
      common /newat/ ia1, iaref(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /files/ imas, nout, ixp, idb, icon
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      do i=1,ia1
        iat = iaref(i)
        write(nout,2) i, atom(iat), 
     *   ( atom( iaref( icmat(i,k) ) ), k=1,nba(i) ) 
      enddo
2     format(i7,1x,a8,' -',10(1x,a8))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c ring perception and test of rings for planarity
c we operate on atom array with all open acyclic atoms already removed
      subroutine FindRings
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /choices/ nbaf(nat,nbond)   
      logical lspe
      dimension idum(nbond)
      common /ring1/ nrmin, nrmax, lspe
      common /newat/ ia1, iaref(nat)
      common /files/ imas, nout, ixp, idb, icon
      dimension id(maxmem), ird(maxr)
      common /visited/ ivis(nat)
      logical aplan
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
c
      nor=0
      do i=1,maxrings
        noar(i)=0
        planar(i)=.false.
        ipldum(i)=0
        do j=1,maxmem
           nra(i,j)=0
        enddo 
      enddo
      call izero(ird,maxmem+1)
      call izero(ivis,nat)
      iatot = ia1
c
      write(nout,'(/80a)') ('-',i=1,80)
      write(   *,345) nrmin, nrmax
      write(nout,345) nrmin, nrmax
      write(*,*) ' '
345   format(/' AUTOMATED SEARCH FOR ',i3,' -',i3,' - MEMBERED RINGS')
c     if there are no atoms in "reduced" structure - skip all this
      if( ia1.le.2 ) then
        write(   *,346) 
        write(nout,346) 
        return      
      endif
346   format(/' THERE ARE NO ATOMS IN THE REDUCED CONNECTIVITY TABLE',
     * //,' SEARCH FOR RINGS WILL NOT BE PERFORMED...')
c
c.....at this stage we have 'ia1' atoms in the atom list
c.....the connectivity matrix is stored in arrays 'icmat' and 'nba'
c.....1st atom in the atom list - 1st atom in test list
c
c McM try to start with atom w/more than 2 connections to 
c     guarantee short rings 
      do i=1,ia1
cDBG         write(nout,*) ' i, nba(i) = ',i,nba(i)
         if (nba(i).gt.2) then 
cDBG         write(nout,*) ' swapping 1 and ',i   
c           swap index 1 and i
            iaref1=iaref(1)
            iaref(1)=iaref(i)
            iaref(i)=iaref1
cDBG            write(nout,*) 'iaref =',(iaref(ij),ij=1,ia1)
c           save 1
	    jdum = nba(1)
	    do j=1,jdum
	      idum(j) = icmat(1,j)
	    enddo
c           copy i-th to 1
	    nba(1) = nba(i)
	    do j=1,nba(i)
	      icmat(1,j) = icmat(i,j)
	    enddo
c           copy saved 1 to i-th	  
	    nba(i) = jdum
	    do j=1,jdum
	      icmat(i,j) = idum(j)
	    enddo
c           now go thru all icmat and exchange 1 and i
	    do j=1,ia1	    
	      do k=1,nba(j)
	        jmat = icmat(j,k)
		if( jmat.eq.1 ) icmat(j,k) = i
		if( jmat.eq.i ) icmat(j,k) = 1
	      enddo
	    enddo
          goto 350
         endif
      enddo
350   continue
c
c     nbaf arrays store info for each connection : 1 - not explored, 0 - explored      
c     initialize nbaf arrays of all atoms in the list
      do i=1,ia1
        do j=1,nba(i)
          nbaf(i,j)=1
        enddo
      enddo   
c     1st atom and its first neighbour
      ird(1) = 1
      ird(2) = icmat(1,1)
c     label as "visited" 
      ivis(1) = 1
      ivis(2) = ird(2)
c     close path to and from 1st atom to its inei-th neigbor
      nbaf(1,1) = 0
c      iat = icmat(1,1)
c     also close it from the "other" side 
c      do j=1,nba(iat)
c        if( icmat(iat,j).eq.1 ) nbaf(iat,j) = 0
c      enddo
c          
100   continue
      ii = 2
111   continue
      il =  ichoice(ird(ii-1),ird(ii))
c       there is path where we can go 
        if( il.ne.0 ) then
          ii = ii + 1
          ird(ii) = il
c         label as visited
          ivis( il ) = 1
cDBG          write(nout,*) '  '
cDBG          write(nout,*) '  ii = ',ii
cDBG          write(nout,*) '     ird = ',(k,k=1,ii)
cDBG          write(nout,*) '       # = ',(ird(k),k=1,ii)
cDBG          write(nout,*) '           ',(atom(iaref(ird(k))),k=1,ii)
c         test for repeated index in array 'ird'
          do j=ii-1,1,-1
cDBG            write(nout,*) ' j, ird(j), il = ',j, ird(j), il
            if( il.eq.ird(j) ) then
cDBG              write(nout,*) ' calling CheckRing',ii,j
              call CheckRing(ii,ird,j,ii)
              goto 111
            endif  
          enddo 
c       go back to the last available choice point
        else
          call Go2LCP(ii,ird,iflag)
c         no more LCPs left - check if some of the atoms were not visited
          if( iflag.eq.0 ) then
             call CheckVisited(iflag,ird)
c            no more "unvisited" atom 
             if( iflag .eq. 0 ) goto 300
c            found "unvisited" atom 
             goto 100
          endif
        endif    
c       add another atom
        goto 111 
c
c.....that's it - no more atoms left in the list
c
300   continue
      write(nout,'(/1x,79a)') ('-',j1=1,79)
      write(nout,301) nor    
      write(   *,301) nor    
301   format(/' TOTAL RINGS FOUND (nor)  = ',i6)
c
c      stop
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c checks if any atoms in the list were not visited
c returns iflag = 0 - ALL atoms were visited
c             <>  0 - found atom which was not visited, then it fills 
c                     in the first 2 members of array 'ird' 
      subroutine CheckVisited(iflag,ird)
      include 'lsdb.inc'
      common /visited/ ivis(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /newat/ ia1, iaref(nat)
      common /choices/ nbaf(nat,nbond)   
      dimension ird(maxr)
c
      iflag = 0
      call izero(ird,maxmem+1)
      
      do 900 i=1,ia1
       if( ivis(i).eq.1 ) goto 900
c      ok, select "unvisited" atom
       ird(1) = i
       do j=1,nba(i)
c        select its "unvisited" neighbor   
         if( nbaf(i,j).eq.1 ) then
           ird(2) = icmat(i,j)
c          close path           
           nbaf(i,j) = 0  
c          misc
           iflag = 1
           ivis(i) = 1
           ivis( ird(2) ) = 1
           return
         endif
         call Err('CheckVisited','Internal error')         
       enddo
900   continue
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c podozrenie na kol'tzo
      subroutine CheckRing(ii,ird,i1,j1)
      include 'lsdb.inc'
      dimension ird(maxr)
      common /files/ imas, nout, ixp, idb, icon
      common /newat/ ia1, iaref(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      logical lspe
      common /ring1/ nrmin, nrmax, lspe
      logical aplan
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      dimension id(maxmem)
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      ir = 0
      do k=i1,j1-1
        ir = ir + 1
        if( ir.gt.nrmax ) return
        id(ir) = iaref(ird(k))
      enddo
      if( ir.lt.nrmin ) return
      if( ird(i1+1).eq.ird(j1-1) ) return
      
      if( ichk_ring(ir,id).eq.1 ) return
c     ok, this is a new ring
      nor = nor + 1           
      if(nor.gt.maxrings) 
     * call Err('FindRings','Too many rings (maxrings)')
      noar(nor) = ir
      do j=1,ir
        nra(nor,j) = id(j)
      enddo  
c     printout some info    
      write(nout,939) nor, noar(nor)
      write(nout,941) (atom(nra(nor,j)),j=1,noar(nor))           
      write(nout,942) (nra(nor,j),j=1,noar(nor))           
939   format(//' **********  FOUND RING # ',i4,' ---> ',
     *       i3,' - MEMBERED  **********'/) 
940   format(3x,99i4) 
941   format(99(2x,a))       
942   format(99(i6,4x) )       
c.....let's assume that the rings is planar (3 is always planar)
      planar(nor)=.true.         
c.....PLANARITY TEST !!!!! for >3 membered rings
      if(noar(nor).ge.4) then
        call Planarity(1,aplan,ir,id,2)
        planar(nor)=aplan
      endif  
      
c900   continue      
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c goes back to the Last Choice Point (LCP) 
c iflag = 0 - no last choice point found
c iflag = i - last LCP we are going to...
      subroutine Go2LCP(ii,ird,iflag)
      include 'lsdb.inc'
      dimension ird(maxr)
      common /files/ imas, nout, ixp, idb, icon
      common /newat/ ia1, iaref(nat)
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /choices/ nbaf(nat,nbond)   
c
      iflag = 1
cDBG      write(nout,'(6x,a)')'GotoLCP :: '
c 
      do 500 j=ii-1,1,-1
c       j-th atom in array ird  
        jat = ird(j)
        if( nba(jat).le.2 ) goto 500
        do k=1,nba(jat)          
          if( nbaf(jat,k).eq.1 ) goto 200
        enddo
500   continue
c     no more choice points left
      iflag = 0
      return
c     ok, we have available choice point in the list
200   continue
      ii = j
cDBG      write(nout,*) ' ii, ird(ii) = ',ii,ird(ii)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c Planarity test for rings/groups. The kind of tolerance to be used
c (i.e. tol2 [ring] or tol3 [group] is given by parameter 'ipt' :
c   ipt =  2 : tol2 [ring]
c       <> 2 : tol3 [group]    
      subroutine Planarity(iprt,planar,imem,id,ipt)
      include 'lsdb.inc'
      logical planar
      dimension xd(3), a(3,3), eval(3), evec(3,3), e(3), q(3)
c     number of atoms in the ring and their cartesian coordinates
      dimension id(maxmem)
      common /files/ imas, nout, ixp, idb, icon
      character atom(nat)*8, asym(nat)*5, alab*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /tols/ tol1, tol2, tol3
      common /verbosity/ iverbose
      data zero /0.00D00/
c.....choose what kind of tolerance to use      
      if( ipt.eq.2 ) then
        toler = tol2
        alab = 'ring '
      else
        toler = tol3
        alab = 'group'
      endif
c
      if(iprt.eq.1) write(nout,278) imem, alab(:ils(alab)), 
     *                              (atom(id(k)),k=1,imem)
278   format(/' Planarity test for ',i3,'-atom ',a,' : ',//,6(2x,a))      
      planar = .true.
c.....calc centroid (or the centre of mass) of points         
      do k=1,3
         xd(k)=zero
         do i=1,imem
            xd(k)=xd(k)+xc(id(i),k)
         enddo
         xd(k)=xd(k)/imem
      enddo       
c.....clear the matrix  
      call rzero(a,9)     
c.....fill in the matrix
      do k=1,imem
        idk = id(k)
        do 301 i=1,3
          do 301 j=1,3
301     a(i,j) = a(i,j) + (xc(idk,i)-xd(i))*(xc(idk,j)-xd(j)) 
      enddo                   
c.....find eigenvalues and eigenvectors using TQL2 subroutine:
c.....1. convert to Tridiagonal
      call Tred2(3,3,a,eval,e,evec)
c.....2. find eigenvalues and eigenvectors
      call Tql2(3,3,eval,e,evec,ierr)
c      write(nout,'(/5x,a,i5)') 'ierr = ',ierr
c.....
c      write(nout,'(/5x,a)') 'Eigenvalues in ascending order:'
c      write(nout,'(3(5x,e15.5/))') eval
c      write(nout,'(5x,a)') 
c     * 'Eigenvectors (n-th column corresponds to the n-th eigenvalue):'
      d=zero
      do i=1,3
        q(i)=evec(i,1)
c        write(nout,'(5x,3e17.7)') (evec(i,j),j=1,3)
        d = d + q(i)*xd(i)
      enddo
      if( iprt.eq.1 .and. iverbose.ge.2 ) write(nout,389) q, d
389   format(/5x,'P Q R S = ',4f12.4/)  
c.....dist to plane for each atom and sigma for plane
      sigpln=zero
      do i=1,imem
        idi = id(i)
C---    (m).(x) - d 
        dist = q(1)*xc(idi,1) + q(2)*xc(idi,2) + q(3)*xc(idi,3) - d
        if( iprt.eq.1 .and. iverbose.ge.2 ) 
     *    write(nout,390) atom(id(i)), dist
390     format(7x,'dist. from ',a,' to lsq plane =',f11.3,' Ang')
        sigpln = sigpln + dist*dist
      enddo           
      sigpln = sqrt( sigpln / (imem-3) )
      if( sigpln.gt.toler ) planar=.false.

      if( iprt.eq.1 ) then
        write(nout,391) sigpln, toler
        write(nout,392) alab(:ils(alab)), planar
      endif  
391   format(/5x,'Esd atoms from the plane (sigpln) = ',f12.4,
     *           ' ( tol = ',f8.4,' )') 
392   format(/5x,'This ',a,' is planar = .',l1,'.')                               
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c* =====================================================================
c* NIST Guide to Available Math Software.
c* Fullsource for module TRED2 from package CMLIB.
c* Retrieved from ARNO on Fri Dec  3 20:15:42 1999.
c* =====================================================================
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C***BEGIN PROLOGUE  TRED2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
C            matrix using and accumulating orthogonal transformation
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED2,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains the real symmetric input matrix.  Only the
C          lower triangle of the matrix need be supplied.
C
C     On Output
C
C        D contains the diagonal elements of the tridiagonal matrix.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C
C        Z contains the orthogonal transformation matrix
C          produced in the reduction.
C
C        A and Z may coincide.  If distinct, A is unaltered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRED2
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE
C
C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0E0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0E0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0E0
      E(1) = 0.0E0
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0E0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0E0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0E0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0E0
            Z(J,I) = 0.0E0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
      END
c
c* ======================================================================
c* NIST Guide to Available Math Software.
c* Fullsource for module TQL2 from package SLATEC.
c* Retrieved from CAMSUN on Fri Dec  3 19:47:06 1999.
c* ======================================================================
c*DECK TQL2
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
C***BEGIN PROLOGUE  TQL2
C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TQL2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional REAL array, dimensioned
C          E(N).
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1, 2, ..., IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TQL2
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL D(*),E(*),Z(NM,*)
      REAL B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0E0
      B = 0.0E0
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0E0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0E0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * E(I) * R
            S = 1.0E0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
*DECK PYTHAG
      REAL FUNCTION PYTHAG (A, B)
C***BEGIN PROLOGUE  PYTHAG
C***SUBSIDIARY
C***PURPOSE  Compute the complex square root of a complex number without
C            destructive overflow or underflow.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PYTHAG-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Finds sqrt(A**2+B**2) without overflow or destructive underflow
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PYTHAG
      REAL A,B
C
      REAL P,Q,R,S,T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cKNJ
c determine if k-th atom belongs to a ring : 
c cyclic = .true. if atoms belongs to at least 1 ring
c plane  = .true. if the ring to which atom belongs if planar         
c iring = number of the ring to which atom belongs (if only 1 ring)
c or sum of the ring numbers if two rings
c plane2 = .true. if the two rings to which atom belongs are planar
c mem = number of ring members
      subroutine IsItCyclic(k,cyclic,plane,plane2,iring,mem)
      include 'lsdb.inc'
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      logical cyclic, plane, plane2
c
      cyclic = .false.
      plane  = .false.
      plane2 = .false.
      iring  = 0
      mem = 0

c.....loop over all found rings
      ij = 0
      do i=1,nor
        do j=1,noar(i)
          if( nra(i,j).eq.k ) then
            ij = ij + 1
c...........number of the 1st ring to which k-th atom belongs and number of its members
            if( ij.eq.1 ) then 
                          nring = i
                          mem1 = noar(i)
            endif 
c...........numbers of the 2 rings to which k-th atom belongs and the sum of their members
            if( ij.eq.2 ) then
                          nring2 = i
                          nring1 = nring
                       nring3 = nring1 + nring2
                          mem2 = noar(i) + mem1
           endif
          endif
        enddo
      enddo
c
c.....Ok, now figure out what to do:
c.....1. Atom doesn't belong to any of the rings
      if( ij.eq.0 ) return

c.....Atom belongs to at least 1 ring
      cyclic = .true.

c.....2. Atom belongs to more than 2 rings
      if( ij.ge.3 ) return

c.....3. Atom belongs to only 1 ring
      if( ij.eq.1 ) then
c.......ring is planar
        if( planar(nring) ) then
          plane  = .true.
          iring = nring
          mem = mem1
          return
        endif
         return  
      endif 
c.....4. Atom belongs to 2 rings
      if( ij.eq.2 ) then
c......both rings are planar
       if(planar(nring1).and.planar(nring2)) then
          plane2  = .true.
          mem = mem2
          iring = nring3
         return
       endif 
c......one of the 2 rings is planar
c.......--> treatment as if it belonegd to one ring only
       if(planar(nring1).and. .not. planar(nring2)) then
          plane  = .true.
          mem = mem1
          iring = nring1
          return
       endif 
       if(planar(nring2).and. .not. planar(nring1)) then
          plane  = .true.
          mem = mem2 - mem1
          iring = nring2
         return
       endif
         return 
      endif
c.......ring is NOT planar
c                    
      return
      end
cKNJ
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c find neigbors for each atom and stores them in common block /neighbors/
c
      subroutine FindAllNeighbors
      Include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /files/ imas, nout, ixp, idb, icon 
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      common /iprint/ iprint  
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
c.....neighbor info from FindNei
      dimension iatn(maxnei), rnd(maxnei) 
      dimension ndt(maxnei), natnum(maxnei), ndtnum(maxnei,maxnei)
      common /neig/ inei, iatn, rnd, ity, ndt, natnum, ndtnum
c
      logical la
      data zero /0.00D00/      
c
      iprint = 0
      write(nout,'(/1x,79a)') ('-',j1=1,79)
      write(   *,734) 
      write(nout,734) 
734   format(/' IDENTIFYING COVALENTLY-BONDED NEIGHBORS FOR EACH ATOM') 
c
c     -------------------
c     loop over all atoms
c     -------------------    
c     
      do 900 i=1,natoms
        call FindNei(i,zero) 
c       fill in common block /neighbors/
        nnei(i) = inei
        do j=1,inei
          jnei(i,j) = iatn(j)
          rnei(i,j) = rnd(j)
        enddo
        jtnei(i) = ity
        do j=1,ity
          janei(i,j) = natnum(j)
          jnnei(i,j) = ndt(j)
          do k=1,ndt(j)
            jsnei(i,j,k) = ndtnum(j,k)
          enddo
        enddo
900   continue
c
c     now for all atoms and all atom neighbor types,
c     sort neighbors within the given type by 
c     1) hybridization/ascending and 2) distance/ascending
c     remember, that we already sorted neighbor types by atomic
c     number 
c
c     -------------------
c     loop over all atoms
c     -------------------    
c     
      write(nout,800) 
800   format(/' Printout of neighbors of all atoms sorted by ',/,
     *  4x,'1) atomic number (descending)',/,
     *  4x,'2) hybridization ( ascending)',/,
     *  4x,'3) distance      ( ascending)')  
      do 902 k=1,natoms
      write(nout,'(/a,i6,3a)')  ' Atom # ',k,' = ',atom(k), ' :'
      do i=jtnei(k),1,-1
c       skip sorting when only 1 neighbor of given type
        if( jnnei(k,i).le.1 ) goto 901
10      la=.false.
        do 5 j=1,jnnei(k,i)-1
          jatI  = jsnei(k,i,j)
          jatI1 = jsnei(k,i,j+1)
c         first check hybridization
          if( nnei(jatI1).lt.nnei(jatI) ) then
            idum = jsnei(k,i,j)
            jsnei(k,i,j) = jsnei(k,i,j+1)
            jsnei(k,i,j+1) = idum
            la = .true.
c         if hybrid. is the same - check dist
          else if( nnei(jatI1).eq.nnei(jatI) ) then
c          get distances for jatI and jatI1 
           do l=1,nnei(k)
             if( jnei(k,l).eq.jatI  ) rI  = rnei(k,l)
             if( jnei(k,l).eq.jatI1 ) rI1 = rnei(k,l)
           enddo  
c          check distances
           if( rI1 .lt. rI ) then
             idum = jsnei(k,i,j)
             jsnei(k,i,j) = jsnei(k,i,j+1)
             jsnei(k,i,j+1) = idum
             la = .true.
           endif            
          endif
5       continue
        if( la ) goto 10
901   continue
c     printout sorted array
      do j=1,jnnei(k,i)
        jat = jsnei(k,i,j)
c       atomic number
        ian = itype(jat)
c       hybridization                
        isp = nnei(jat)
c       get distance to k-th atom
        do l=1,nnei(k)
          if( jnei(k,l).eq.jat ) rd = rnei(k,l)
        enddo
c       print it
        write(nout,801) ian, isp, rd, jat, atom(jat)
801     format(5x,2i4,f11.3,i8,3x,a)                
      enddo
c        
      enddo        
902   continue
c    
c     printout
c
c      do i=1,natoms
c        call PrtAtomNei(i,nout)
c      enddo
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine LocalDef
      include 'lsdb.inc'
c
      character line*200, atom(nat)*8, asym(nat)*5, key*68
      dimension xf(nat,3), xc(nat,3), itype(nat)
      dimension xj(3)
      character dlab*6, xdold*60
      common /atoms/ ia, atom, xf, xc, itype, asym
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /files/ imas, nout, ixp, idb, icon
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      common /kappaset/ itkset, kset(nat)      
      common /kappas/ nkap,rkappa(nat,2),nkapsf(nat)
      common /ftype/ iform
      common /ifra/ ifra
      common /xdold/ xdold
      logical lsx
      common /lsx/ lsx
      logical lkeys
      common /lkeys/ lkeys(7)
      
!! to comment
      character dstring*30
      character atom1*10
      character atom2*10
c.....
      ndum = 0
c.....check: create PLATON file with cartesian coordinates
      iz=52
      open(iz,file='lsdb.ang',form='formatted',status='unknown')
      write(iz,'(a)') 'ANGSTROM'
53    format('ATOM ',a,3x,3f14.4)      
c
c....if read XD files - create copy from old to new
c
      if(iform.eq.2) then
        iold=39
        open(iold,file=xdold,form='formatted',status='old')
123     read(iold,'(a)') line
c         write(*,999) 'line=[',line,']'
         if(line(1:4).ne.'ATOM') then
c           write(*,*) 'ils=',ils(line)
           write(imas,'(a)') line(1:ils(line))
           goto 123 
         endif            
      endif 
c      
      write(   *,  *) '   '
      write(   *,678)
      write(imas,678)  
678   format('ATOM     ATOM0    AX1 ATOM1    ATOM2   AX2 R/L TP TBL KAP 
     *LMX SITESYM  CHEMCON')     
      do 900 i=1,ia
         write(nout,'(/1x,79a)') ('-',j1=1,79)
         write(nout,800) i, atom(i), occ(i)
800      format(/' Working with atom ',i7,' : [',a,']',5x,'occ =',f11.4)
        call DefLoc(i)
c     
c         write(*,*) ' i = ',i,'   asym(25) = ',asym(25)
c........write atom to PLATON file
         write(iz,53) atom(i), (xc(i,j),j=1,3)
900   continue
c.....write dummy atoms
      do i=1,ndum
        write(imas,'(a6,2x,3f12.6)') dlab(i), (xdum(i,k),k=1,3)
c.......cartesian coordinates to PLATON file
        if(ifra.eq.1) then
          call rabs(0,xdum(i,1),xdum(i,2),xdum(i,3),xj(1),xj(2),xj(3))
          write(iz,679) i, xj
        else
          write(iz,679) i, (xdum(i,k),k=1,3)
        endif  
      enddo
679   format('ATOM Li(',i1,')',2x,3f14.4)      
c.....finished ATOM table      
      write(imas,999) 'END ATOM'
c.....finished PLATON file
      close(iz) 
c
c.....iform=2: xd files
c  
      if( iform.eq.2 ) then    
345     read(iold,'(a)') line
        call UpCase(line)
        if(line(1:4).ne.'END ') goto 345
346     read(iold,'(a)') line
        if(line(1:4).ne.'KEY ') then
           write(imas,'(a)') line(1:ils(line))
           goto 346
        endif
      endif  
c
c......shelx --> xd.mas 
c
      if( iform.eq.1 ) then
        call WrXDmas1(imas,itkset)
      endif
c
c.....write out the refined parameters based on symmetry
c
c      write(imas,999) '!'
      write(imas,999) 'KEY     XYZ --U2-- ----U3---- ------U4------- M- 
     *-D- --Q-- ---O--- ----H----'
      do 901 i=1,ia
        call GetKey(i,key)
        write(imas,902) atom(i),key
901   continue      
902   format(a,a)
      do 4 i=1,itkset
c       if just create XD files w/o DB 
c       or kappas were not requested
        if( lsx .or. .not.lkeys(7) ) then
          write(imas,999) 'KAPPA   000000'    
          goto 4    
        endif
        do k=1,ia
c........non-H atoms        
          if( kset(k).eq.i .and. itype(k).ne.1) then
            write(imas,999) 'KAPPA   110000'        
            goto 4
          endif 
        enddo   
c.......H-atoms        
        write(imas,999) 'KAPPA   100000'        
4     continue
c
c.....SHELX or 4TOPOND
      if(iform.eq.1.or.iform.eq.3) then
        write(imas,999) 'EXTCN   0000000'
        write(imas,999) 'OVTHP   0'
        write(imas,999) 'SCALE   1'
        write(imas,999) 'END KEY'
c.....XD file
      else if(iform.eq.2) then
347     read(iold,'(a)') line
        call UpCase(line)
        if(line(1:6).ne.'EXTCN ') goto 347
        write(imas,'(a)') line(1:ils(line))
348     read(iold,'(a)',end=349) line
          write(imas,'(a)') line(1:ils(line))
        goto 348
349     close(iold)
      endif              
c
c.....SHELX      
c
      if( iform.eq.1 ) then
        write(imas,999) '!',('-',i=1,78)
        write(imas,999) '!'
        write(imas,999) '   END XDLSM'
        write(imas,999) '!'
        write(imas,999) ('!',i=1,79)
        call WrXDmas2(imas)
      endif
c      
999   format(100a)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
      subroutine GetKey(k,key)
      include 'lsdb.inc'
      logical ichemcon
      character atom(nat)*8, asym(nat)*5, mkey*30, skey*38, key*68     
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /chemcon/ ichemcon(nat)
      logical lsx
      common /lsx/ lsx
      logical lkeys
      common /lkeys/ lkeys(7)
c.....keys for structural parameters
      skey = '000 000000 0000000000 000000000000000 '
c.....keys for multipoles
      mkey = '00 000 00000 0000000 000000000'
c....."total" key
      key =  skey//mkey
c
c     figure our what to do with structural keys
c      
      if( .not.lkeys(1) ) goto 10
      if( lkeys(2) .or. (lkeys(3).and.iaf(k).ne.-1) ) then
        skey =                 '111 111111 0000000000 000000000000000 '
c       hydrogen
        if(itype(k).eq.1) skey='111 100000 0000000000 000000000000000 '
      endif
c            
10    continue
      key = skey//mkey
c
c     conversion from shelx file - no multipoles or
c     did not ask for multipoles
c
      if( lsx .or. .not.lkeys(4) ) return
c
c     MULTIPOLE PARAMETERS 
c      
      if( lkeys(5) .or. (lkeys(6).and.iaf(k).ne.-1) ) then
c                                 '00 000 00000 0000000 000000000'
        if(asym(k).eq.'NO')   mkey='10 111 11111 1111111 111111111'      
        if(asym(k).eq.'C')    mkey='10 000 11111 0000000 111111111'      
c
        if(asym(k).eq.'2')    mkey='10 000 10011 1001100 100110011'      
        if(asym(k).eq.'m')    mkey='10 110 10011 0110011 100110011'      
        if(asym(k).eq.'2/m')  mkey='10 000 10011 0000000 100110011'      
c
        if(asym(k).eq.'222')  mkey='10 001 10010 1000100 100100010'      
        if(asym(k).eq.'mm2')  mkey='10 001 10010 1001000 100100010'     
        if(asym(k).eq.'mmm')  mkey='10 000 10010 0000000 100100010'     
c
        if(asym(k).eq.'4')    mkey='10 001 10000 1000000 100000011'     
        if(asym(k).eq.'-42m') mkey='10 001 10000 0001000 100000010'     
c       
        if(asym(k).eq.'3')    mkey='10 001 10000 1000011 100001100'     
        if(asym(k).eq.'3m')   mkey='10 001 10000 1000010 100001000'     
c      
        if(asym(k).eq.'cyl')  mkey='10 001 10000 1000000 100000000'
c
      endif
      key = skey//mkey
c.....      
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  
      subroutine WrXDmas1(imas,nkap)    
      common /shelx1/ fvar(200), wght(6), ifv
90    format(100a)
91    format( a, 20i4,' \', 10(/,20i4,' \') )
      write(imas,90) '!',('-',i=1,78)
      write(imas,91) 'KEEP  kappa ',(i,i=1,nkap)
      write(imas,90) 'KEEP  charge group1'
      write(imas,90) '!KEEP  rigid group1'
      write(imas,90) '!RESET bond C(1) H(1) 1.09 ...'
      write(imas,92) 'WEIGHT  ',wght
92    format(a,6f11.5)      
      write(imas,90) 'SKIP  obs 0. 1.d10 *sigobs 3. 1.d06 sinthl 0. 2.'
      write(imas,90) 'PRINT sinthl .0 2. obs 0. 15. delta 0. 10. *del% 8
     *0 100 extcn 80. 100. *abssc'
      write(imas,90) '!EXTCN  *iso aniso *type_1 type_2 type_3 distr_g *
     *distr_l msc_0  msc_1'
      write(imas,90) '!DMSDA  1.1  1.8'
      write(imas,90) '!FOUR  fmod1 4 2 0 0  fmod2 -1 2 0 0'
      write(imas,90) '!CON  num1 par1/iat1 num2 par2/iat2 ... = num0'
      write(imas,90) '!',('-',i=1,78)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     
c**** Module XDProp
c     
      subroutine WrXDmas2(nm) 
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
c
 99   format(79('!'))
 399  format(100a)
 398  format(a,i4,a)
      write(nm,399) '!'
      write(nm,399) '   MODULE *XDPROP'
      write(nm,399) '!'
      write(nm,399)'MODEL  iam *multipole'
      write(nm,399)
     1     '!APPLY  symm 1 translations 0 0 0 ato(1) ato(2) ...'
      write(nm,399)'SELECT  *local  numdx  check  esd  nocore'
      write(nm,399)
     1     'SELECT  cpcut 0.0001 lmax 4  nstep 20  rcut 4.0 '
      write(nm,399)'SELECT  scale 0.05  dx 0.001  ds 0.005 '
      write(nm,399)
     1     'SELECT  rad1 0.1  rad2 200.  rad3 10.  zone1 1  zone2 1'
      write(nm,399)'SELECT  rEcrit 0.00001  rNcrit 0.00001'
      write(nm,399)
     *'QUADINT   iqt 2  Nrad 50  Nang 194  Becke  *Stock  Print'
c
      write(nm,399)'!GROUP  ato(1) ato(2) ...'
      write(nm,399)'MOLMOM  *cmass center ucell ccharge'
      write(nm,399)'!D-POP'
      write(nm,399)'!EXPORT *orient *min16 lmax 4 nmol 1 natmol ...'
      write(nm,399)'PROPERTY  *rho gradrho d2rho nucpot elpot core', 
     1 ' valence defden sigrho siglap esp ef efg'
      write(nm,399)'!QFIT  grid 11 length 7.0 width 1.0 constrain false'
      write(nm,399)'!CONSTRAIN  ato(1) ato(2) ...'
      write(nm,399)'!'
      write(nm,399)'!STOCKMOM  defden lmin 0 lmax 4 *cmass center ucell  
     *debug'
      write(nm,399)
     * '!STOCKMOM *orient *min16 lmax 4 nmol 1 natmol ...'
      write(nm,399)'!STOCKMOM  atoms *all select ato(1) ato(2) ...'  
      write(nm,399)'!'
      write(nm,399)'!ATATPOT  *EXREP spack  wilcox lj *saptP  saptA '
      write(nm,399)'!ATATPOT  *DISPR spack  wilcox lj *saptP  saptA '
      write(nm,399)'!ATATPOT  *INDUC                  *saptP  saptA'
      write(nm,399)'!'
      write(nm,399)'!HBONDS   H(1) O(1) 1.75   H(2) O(2) 1.82'
      write(nm,399)'!'
      write(nm,399)'!HPOLAR   H(1) H(3) H(4)'
      write(nm,399)'!'
      write(nm,399)'!INTEREN  frag 1  nat1 -nat2  *neutral'
      write(nm,399)'!INTEREN  frag 2  nat3 -nat4  *neutral'
      write(nm,399)
     *'!INTEREN  EP aMM  mMM *EPMM  rCrit1 4.5 rCrit2 20. debug'
      write(nm,399)'!'
      write(nm,398)'!LATEN  frag1 1 ',-ia,'  *neutral'
      write(nm,399)'!LATEN  ncells 0 20  nappl 0  symms 1 2 3'
      write(nm,399)
     *'!LATEN  EP  aMM  mMM  *EPMM  rCrit1 4.5 rCrit2 20. lapf debug'
      write(nm,399)'!'
      write(nm,399)'!POINT  x y z'
      write(nm,399)'!LINE  ato(1) ato(2)  npts 50'
      write(nm,399)'!LINE  points x1 y1 z1 x2 y2 z2 npts 50'

      write(nm,399)'!CUBE  centre x y z   npts 30  stepsize 0.1'
      write(nm,399)'!CUBE  ato(1) ato(2)  npts 20  stepsize 0.1'
      write(nm,399)
     1 '!MAP  atoms ato(1) ato(2) ato(3) npts 50 stepsize 0.1'
      write(nm,399)
     1 '!MAP  bvect1 x1 y1 z1 bvect2 x2 y2 z2 cen x0 y0 z0',
     2 ' npts 50 stepsize 0.1'
      write(nm,399)'!CPSEARCH  bond ato(1) ato(2)'
      write(nm,399)'!CPSEARCH  bond rmin  1.2 rmax  1.6'
      write(nm,399)'!CPSEARCH  ring ato(1) ato(2) ...'
      write(nm,399)
     1     '!CPSEARCH  shell ato(1) rmin 0.3 rmax 0.5 nrad 5 ', 
     2     'nang 11 11 cutoff 16.0'
      write(nm,399)
     1'!CPSEARCH BUBBLE ato(1) rmin 0.3 rmax 0.5 curv -3 ncps 3'
      write(nm,399)'!CPSEARCH  point x y z'
      write(nm,399)'!CPSEARCH  start file.cps'
      write(nm,399)'!BPATH  ato(1) ato(2) algrithm 2'
      write(nm,399) '!'
      write(nm,399)'   END XDPROP'
      write(nm,399)
     1 '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     2 '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(nm,399)'   MODULE XDGEOM'
      write(nm,399)
     1     'SELECT   rmin 0.8 rmax 1.8  tor *ato *bon *ang loc non'
      write(nm,399)'   END  XDGEOM'
      write(nm,99)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c
c     ACTUALLY DEFINE THE LOCAL COORDINATE SYSTEM FOR k-th ATOM
c
c
      subroutine DefLoc( k )
      include 'lsdb.inc'
c
      character*90, dimension(nat,1) :: MOPRO
      common /mopro/ MOPRO
      character*90, dimension(nat,1) :: MCON
      common /mcon/ MCON
      character dstring*30
      dimension xca1(3), xca2(3), xm(3)
      character ax1, ax2, at1*8, at2*8, at*8
      character atom(nat)*8, asym(nat)*5
      character atom1*10, atom2*10
      character ato1*10, ato2*10, ato*10, ato3*10, ato4*10
      character sym*4
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      dimension id(maxmem)
      common /isn/ isn(maxnei)
      common /simple/ iSimple
      logical resi11
      character*17, dimension(nat,1):: RESI1
      character*10, dimension(nat,1):: RESI2
      common /resi1/ RESI1, RESI2, resi11
c.....neighbor info from FindNei
      dimension iatn(maxnei), rnd(maxnei) 
      dimension ndt(maxnei), natnum(maxnei), ndtnum(maxnei,maxnei)
      common /neig/ inei, iatn, rnd, ity, ndt, natnum, ndtnum
      common /files/ imas, nout, ixp, idb, icon
      logical planar, SameAtoms
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      logical lspe
      common /ring1/ nrmin, nrmax, lspe
      character dlab*6
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      logical lfH
      common /lfH/ lfH
      logical lreaddb
      logical cyclic, plane, plane2
      logical lPG ! true when atomic group is 'planar enough'
      common /lPG/ lPG
c011602: for kappa sets
      logical chemcon, ichemcon
      common /kappaset/ itkset, kset(nat)      
      common /iprint/ iprint    
      dimension isame(maxnei,2)
      common /chemcon/ ichemcon(nat)
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
c.....exceptions
      character group(nat)*40
      common /except/ iexp(nat), group
      logical ldb
      character axes(3)*1
      common /ldb/ ldb, ntot, nfound, nmissed
      common /mpoles1/ ifind, Pv1, sPv1, rkap1(2), Plm1(0:4,-4:4)
      common /kappas/ nkap,rkappa(nat,2),nkapsf(nat)
      common /mpoles/ Pv(nat), sPv(nat), Plm(nat,0:4,-4:4), iar(nat,11)
      common /forH/ isp(nat), ihp(nat), ihnei(nat)
      common /forRing/ iring, imemb
      common /ftype/ iform
      common /ifra/ ifra
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      logical lhyd, lheav
      common /lchem/ lhyd, lheav
      logical lsx
      common /lsx/ lsx
      logical lkdis, lndis      
      logical led, ldd
      common /disorder/ led, ldd
      common /verbosity/ iverbose
      character line*80
      character line2*80
      character line3*80
      data axes /'X','Y','Z'/
      data zero / 0.00D00 /, one / 1.00D00 /
c
      ax1 = ' '
      ax2 = ' '
      at1 = ' '
      at2 = ' '
      sym = ' '
      iring = 0
      imemb = 0
c.....misc
      do j=1,3
        xca1(j) = zero
        xca2(j) = zero
      enddo
      iprint = 1
c      write(nout,*) ' k = ',k
c      call FindNei( k, zero )      
c      call SaveNei
c      write(nout,*) ' k = ',k
      call PrtAtomNei(k,nout)
c.....get sorted-in-a-special-way list of neighbors
      call SortNei(k)     
c011602      
c      knei = inei  
c.....check if k-th atom or any of its neighbors are disordered
      lkdis = .false.
      if( iaf(k).eq.-1 ) lkdis = .true.  ! k-th atom itself is disordered
      lndis = .false.
      do i=1,nnei(k)
        iat = jnei(k,i)
        if( iaf(iat).eq.-1 ) lndis = .true. ! one of its neighbors is disordered
      enddo
c.....defaults for lPG   
      lPG = .true.
      if( nnei(k).ge.4 ) lPG = .false.
c.....seq. numbers of first two neighbors in /atoms/
      if( nnei(k).ge.1 ) i1st = jnei(k,1)
      if( nnei(k).ge.2 ) i2nd = jnei(k,2)
c.....how many attached H-atoms
      ih = 0
      if( janei(k,1).eq.1 ) ih = jnnei(k,1)
      write(nout,'(/a,i4)') ' Number of covalently-bonded H-atoms = ',ih
c
c ===========================================
c  Just want to create XD files w/o DB 
c ===========================================
c
      if( lsx ) then
        if( nnei(k).ge.2 ) goto 5000
        goto 4000
      endif
c 
c ===============================================================
c  NO covalently-bonded neighbors were found 
c  So we have to find any neighbors to set local coodinate system
c  Then we don't care about sorting them in any special order
c ===============================================================
c
      if( nnei(k).ne.0 ) goto 5000
c
c     try to find neighbors within 'rdum' Angstroms
c
c     start at 5 Angstroms and increase until 20 Angstroms 
      rdum = 5.00D00
4000  continue    
c      write(*,*) ' rdum = ',rdum
      call FindNei( k, rdum )   
      if( inei.lt.2 .and. rdum.lt.20.00D00 ) then
        rdum = rdum + one
        goto 4000
      endif
c.....Ok, found at least 2 neighbors within 'rdum' Angstroms
      if( inei.ge.2 ) then
c        write(nout,'(a,3i4)') ' k check isn(1..2) = ',k,isn(1),isn(2)
        at1 = atom(iatn(1))
        ax1 = 'X'
        at2 = atom(iatn(2))
        ax2 = 'Y'
        asym(k) = 'NO'
        do j=1,3
          xca1(j) = xc(iatn(1),j)
          xca2(j) = xc(iatn(2),j)
        enddo
        goto 71
c.....if no neighbors found within 'rdum' Angstroms :
c.....make dummies at 0,0,0 and 0.5,0.5,0.5               
      else 
        write(   *,384) atom(k), rdum
        write(nout,384) atom(k), rdum
        write(   *,385)
        write(nout,385)
384     format(/' No neighbors of atom ',a,' were found within ',
     *          f10.2,' Angstroms') 
385     format(/' Will make dummies at 0,0,0 and 0.5,0.5,0.5')             
c
        call AddDummyAtom(ifra,0.,0.,0.,idum,xca1(1),xca1(2),xca1(3))
        at1 = dlab(idum)
        ax1 = 'X'
        call AddDummyAtom(ifra,0.5,0.5,0.5,idum,xca2(1),xca2(2),xca2(3))
        at2 = dlab(idum)
        ax2 = 'Y'
        asym(k) = 'NO'
        goto 71
      endif
      call Err('DefLoc','This should not have happened (before 5000)')
5000  continue
c 
c ===========================================================
c  More than 4 covalently-bonded neighbors  or
c  just want to create XD files w/o DB   or
c  one of the neighbors is disordered but atom itself is not
c ===========================================================
c
      if( nnei(k).gt.4 .or. lsx .or. (.not.lkdis.and.lndis) ) then
        at1 = atom(isn(1))
        ax1 = 'X'
        at2 = atom(isn(2))
        ax2 = 'Y'
        asym(k) = 'NO'
        do j=1,3
          xca1(j) = xc(isn(1),j)
          xca2(j) = xc(isn(2),j)
        enddo
        goto 71
      endif  
c
c =================================================
c     DETERMINE IF THE ATOM BELONGS TO A RING
c =================================================
cKNJ
c.....Determine whether this atom is in the ring & if the ring is planar
      call IsItCyclic(k,cyclic,plane,plane2,ir,mem)
c
c ===========================================================
c     Atom belongs to 1 PLANAR ring - mark it or special case
c ===========================================================
c
      if( cyclic .and. plane ) then
        write(nout,'(/a,i5,a)') 
     *   ' Atom belongs to only one planar ring ( # ',ir,' )'
        iring = 1
        imemb = mem
        if( .not.lspe ) goto 827
c       process atoms in 1 planar ring as a special case
c.......check to see if DUM atom for this ring has already been created
        if( ipldum(ir).eq.0 ) then
c.........coordinates of the dummy atom - centroid         
          do l=1,3
            xm(l) = zero
            do i = 1, noar(ir)
              iat = nra(ir,i)
              if(ifra.eq.1) xm(l) = xm(l) + xf(iat,l)
              if(ifra.eq.0) xm(l) = xm(l) + xc(iat,l)
            enddo
            xm(l) = xm(l) / noar(ir)
          enddo  
c.........try to add dummy atom         
          call AddDummyAtom(ifra,xm(1),xm(2),xm(3),idum,xj,yj,zj)
c.........reference to the dummy atom, which defines the center of the ring
          ipldum(ir) = idum 
          kdum = ndum      
        else
          kdum = ipldum(ir)
        endif                
        write(nout,569) dlab( kdum ) 
569     format(/' Dummy atom which defines the ring center = ',a)
c.......OK, dummy atom (centre of the ring) is already defined.
c.......X-axis points towards the centre of the ring:
        asym(k) = 'm'
        at1=dlab( kdum )
        ax1 = 'X'
        if(ifra.eq.1) then
          call rabs(0,xdum(kdum,1),xdum(kdum,2),xdum(kdum,3),
     *                     xca1(1),     xca1(2),     xca1(3) )
        else
          do l=1,3
            xca1(l) = xdum(kdum,l)
          enddo
        endif
c.......Y-axis is directed towards first atom in 'isn' array that belongs
c.......to ir-th ring         
        ax2 = 'Y'
        do j=1,nnei(k)
          jat = isn(j)
          do l=1,noar(ir)
            if( nra(ir,l).eq.jat ) goto 566
          enddo
        enddo
        call Err('DefLoc','Could not find atom to direct Y-axis')
566     continue
        at2 = atom(jat)
        do j=1,3
          xca2(j) = xc(jat,j)
        enddo
        goto 71
c
c ======================================================================
c     Atom belongs to 1 ring but it is NOT PLANAR (enough) or it belongs
c     to more than 2 rings. Anyway it will be treated as NOT_CYCLIC atom
c ======================================================================
c
cKNJ
c
      else if( cyclic .and. .not.plane .and. .not. plane2 ) then
        write(nout,828)         
      endif
c
c ======================================================================
c     Atom belongs to 2 PLANAR rings - mark it or special case
c ======================================================================
c
      if( cyclic .and. plane2 ) then   
        write(nout,'(/a,i5,a)') 
     *   ' Atom belongs to 2 planar rings of tot numb ( #',ir,' )'
       iring = 2
       imemb = mem
        if( .not.lspe ) goto 827
        asym(k) = 'm'
        ax1 = 'X'
        do j=1,nnei(k)
         jat = isn(j)
        call IsItCyclic(jat,cyclic,plane,plane2,ir,mem)
            if( cyclic .and. plane2 ) goto 567
            if( cyclic .and. .not. plane) goto 567 
        enddo
        call Err('DefLoc','Could not find atom to direct X-axis')
567     continue
        at1 = atom(jat)
        do j=1,3
          xca1(j) = xc(jat,j)
        enddo
c
        if( mem.eq.11 .or. mem.eq.12) then
        ax2 = 'Y'
        do j=1,nnei(k)
          jat = isn(j)
        call IsItCyclic(jat,cyclic,plane,plane2,ir,mem)
            if( cyclic .and. plane .and. mem.eq.6) goto 568 
        enddo
        call Err('DefLoc','Could not find atom to direct Y-axis')
        endif
c 
        if( mem.ne.11 .or. mem.ne. 12) goto 827
c
568     continue
        at2 = atom(jat)
        do j=1,3
          xca2(j) = xc(jat,j)
        enddo
        goto 71 
      endif
c
828   format(/' This atom belongs to 1 or 2 rings which are NOT planar ',
     *' or it belongs to >2 rings',/,
     *' Will be treated as if it DOES NOT belong to a ring')
cKNJend
c =====================================
c     ATOMIC GROUP PLANARITY CHECK
c =====================================
c
827   continue
      if( nnei(k).lt.3 ) goto 448
      call CollectNei(k,id,imem)
      idp = 0
      if( iverbose.ge.3 ) idp = 1
      call Planarity(idp,lPG,imem,id,3)
      write(nout,202) imem, lPG 
202   format(/' Atomic group (',i2,' atoms) planarity (lPG) = .',l1,'.')
448   continue
c
c =================================================
c     PROCEED AS USUAL - NOT CYCLIC ATOM
c =================================================
c
c     Let's first check the environment of each of 'inei' neighbors 
c     to figure our what neighbors among 'nnei(k)' are equivalent
c     we do not compare here functional groups and local symmetries...
c
      is = 0  ! number of 'same' pairs
      do 54 i=1,nnei(k)
        do 54 j=i+1,nnei(k)
         iat = jnei(k,i)
         jat = jnei(k,j)
         if( SameAtoms(iat,jat,0,0,iSimple,0) ) then
           is = is + 1
c          isame(is,1..2) - stores is-th pair of "same" atoms
           isame(is,1) = iat
           isame(is,2) = jat
           write(nout,'(/a,2i6)') ' Same atoms : ',iat,jat
         endif
54    continue
      write(nout,'(/a,i4)') ' is = ',is
c....
      goto(10,20,30,40) nnei(k)
c
c ===============================
c  1 neighbor 
c ===============================
c      
10    continue
      if( itype(k).ne.1 ) goto 111
c
c     - H-atom: simple? - not quite! -
c
      ax1='Z'
      at1=atom(i1st)        ! parent atom
      do j=1,3
        xca1(j)=xc(i1st,j)
      enddo
      ax2='Y'
      asym(k) = 'cyl'        
      ihp(k) = itype(i1st)  ! atomic number of parent atom
c     find nei of parent atom 
c      call FindNei( i1st, zero )      
c     sp-hybridization of parent atom
      isp(k) =  nnei(i1st)
c     number of attached H-atoms        
      ihnei(k) = jnnei(i1st,1) 
      write(nout,836) isp(k), ihnei(k)
836   format(/' Hydrogen connected to an atom with ',i2,
     *        ' bonds ( ',i2,' - Hydrogen(s) )')      
c.....now find atom to direct the 2nd axis - must not be H-atom, or
c.....at least not this (k-th) H-atom
c.....still working with environment of parent atom.... 
c
c     only 1 type (H) but at least 2 H-atoms (H2O) :  
      if( jtnei(i1st).eq.1 .and. nnei(i1st).ge.2 ) then
        do i=1,nnei(i1st)
          iat = jnei(i1st,i)
          if( iat.ne.k ) then
            at2 = atom(iat) 
            do j=1,3
              xca2(j)=xc(iat,j)
            enddo
            goto 11
          endif
        enddo 
        call Err('DefLoc','can not happen <1>')
      endif
c     general case - more than two atom types or just 1 H-atom 
c     find neigbors within 10 Ang from parent atom
c     pick up 1st non-H atom but angle should be < 170 deg
c     also make sure that chosen atom is NOT symmetry-generated
      iprint = 0
      call FindNei( i1st, 10. )  
      do i=1,inei
        iat = iatn(i)
        if( itype(iat).ne.1 .and. isymeq(iat).eq.0 ) then
c         check the angle H-parent-nei
          call CalcAng( xc(k,1), xc(k,2), xc(k,3), 
     *                  xc(i1st,1), xc(i1st,2), xc(i1st,3),
     *                  xc(iat,1), xc(iat,2), xc(iat,3), ang  )
          if( ang.le.170.d0 ) then
            call PrtAng(nout,k,i1st,iat,ang)     
            at2=atom(iat)
            do j=1,3
              xca2(j)=xc(iat,j)
            enddo
            goto 11
          endif  
        endif
      enddo
c     ok, failed to find suitable neighbors within 10 Ang
c     make dummy at 0,0,0               
      write(   *,384) atom(i1st), rdum
      write(nout,384) atom(i1st), rdum
      write(   *,386)
      write(nout,386)
      call AddDummyAtom(ifra,0.,0.,0.,idum,xca2(1),xca2(2),xca2(3))
      at2 = dlab(idum)
11    continue
c      
c     finally extend X--H distance to standard/neutron value    
c
      if( lfH ) call ExtendH(i1st,k)
      goto 71
c
386   format(/' Will make dummy at 0,0,0')  
c      
111   continue
c
c    - Other atoms -
c
c  ? - X    group
c  ------------------------------
c  X-axis - along the ?-X bond
c  Y-axis - along ?-* bond, where * is the closest non-H neighbor of X                
c
c.....1st atom is the closest neighbor of k-th atom
      ax1 = 'X'
      at1=atom(i1st)        ! neighbor
      do j=1,3
        xca1(j)=xc(i1st,j)
      enddo
      ax2 = 'Y'
      asym(k) = 'm' 
      i3rd = 0   
c     if there are at least 2 neighbors connected to 'i1st' atom :
      if( nnei(i1st).ge.2 ) then
        do 438 i=1,nnei(i1st)
          i3rd = jnei(i1st,i)
          if( i3rd.ne.k ) then
            call CalcAng( xc(k,1),    xc(k,2),    xc(k,3),
     *                    xc(i1st,1), xc(i1st,2), xc(i1st,3),
     *                    xc(i3rd,1), xc(i3rd,2), xc(i3rd,3), ang ) 
c...........if angle is too small (i.e. close to 0, which means
c...........that i3rd-th atom is close to k-th atom), then we should
c...........choose another one...
            if( ang.le.10.D0 ) goto 438  
c...........bent conformation
            if(ang.le.160.d0) then
              call PrtAng(nout,k,i1st,i3rd,ang)
              at2=atom(i3rd)
              do j=1,3
                xca2(j)=xc(i3rd,j)
              enddo
              goto 71
c...........linear conformation: change axes assignments
            else
              ax1='Z'
              asym(k) = 'cyl'   
            endif  
          endif
438     continue
      endif
c.....ok, this is either a linear conformation or 'i1st' is only bonded
c.....to k-th atom; in any case we have to find another atom in 'bend' 
c.....position to direct the 2nd axis        
      rdum = 4.00D00
12    continue
      call FindNei( i1st, rdum ) 
      do 439 i=1,inei
        ij = iatn(i) 
        if( ij.ne.k .and. (i3rd.ne.0.and.i3rd.ne.ij) ) then
c.........check the angle k-i1st-ij
          call CalcAng( xc(k,1),    xc(k,2),    xc(k,3),
     *                  xc(i1st,1), xc(i1st,2), xc(i1st,3),
     *                  xc(ij,1),   xc(ij,2),   xc(ij,3),ang )      
c.........if angle is too small (i.e. close to 0), then we should
c.........choose another atom
          if( ang.le.10.D0 ) goto 439  
c.........ok, this angle is acceptable
          if( ang.le.160.d0 ) then
            call PrtAng(nout,k,i1st,ij,ang)
            at2=atom(ij)
            do j=1,3
              xca2(j)=xc(ij,j)
            enddo
            goto 71
          endif  
        endif
439    continue
c.....if no neighbors found within ~10 Angs - make dummy at 0,0,0               
      if(rdum.gt.10.0) then
        write(   *,384) atom(i1st), rdum
        write(nout,384) atom(i1st), rdum
        write(   *,386)
        write(nout,386)
        call AddDummyAtom(ifra,0.,0.,0.,idum,xca2(1),xca2(2),xca2(3))
        at2 = dlab(idum)
        goto 71
      endif
      rdum = rdum + 2.D0 ! if no neigbours within rdum are found - increase rdum
      goto 12
c
c ==============================
c  2 neighbors - sp C, O, N etc
c ==============================
c      
20    continue
c
c  i1st - k - i2nd   group 
c  
c.....let's first check the angle i1st-k-i2nd
      call CalcAng( xc(i1st,1), xc(i1st,2), xc(i1st,3),
     *              xc(k,1),    xc(k,2),    xc(k,3),
     *              xc(i2nd,1), xc(i2nd,2), xc(i2nd,3),dang )  
      if( dang.ge.170.0d0 ) goto 21
c
c     bent conformation: "m" or "mm2"
c  
      write(nout,'(/a)') ' BENT conformation !'      
c
c     i1st and i2nd are the same (including H-atoms): mm2
c
c      ndum=0
      if(is.eq.1) then
        call MidPoint( xc(i1st,1),xc(i1st,2),xc(i1st,3),
     *                 xc(i2nd,1),xc(i2nd,2),xc(i2nd,3),
     *                      xm(1),     xm(2),     xm(3) ) 
        call AddDummyAtom( 0, xm(1),xm(2),xm(3), idum, 
     *                     xca1(1),xca1(2),xca1(3) )     
        at1=dlab(idum)
        ax1='Z'
        at2=atom( i1st )
        ax2='Y'
        asym(k) = 'mm2'
        do l=1,3
          xca2(l)=xc(i1st,l)
        enddo
        goto 71
c
c     X and Y are different: m
c
      else
        ax1 = 'X'
        ax2 = 'Y'
        asym(k) = 'm'
c.......special case i1st = H : X-axis is along k-i2nd rather than k-i1st(H)
        if( ih.eq.1 ) then
          i1 = i2nd
          i2 = i1st
c.......general case 
        else   
          i1 = i1st
          i2 = i2nd
        endif
        at1 = atom(i1)
        at2 = atom(i2)
        do j=1,3
          xca1(j) = xc(i1,j)
          xca2(j) = xc(i2,j)
        enddo
        goto 71
      endif  
c
c     linear conformation:  ONLY bond-directed multipoles
c
c     Z-axis along the closest non-H atom, or along the closest
c     H-atom if both neighbors are H-atoms    
c     X-axis goes towards dummy at 0,0,0 if there are no other choices
c
21    continue
      write(nout,'(/a)') ' LINEAR conformation !'      
      ax1 = 'Z'
      ax2 = 'X'
      asym(k) = 'cyl'
      in1 = i1st
      if( itype(i2nd).ne.1 ) in1 = i2nd
      at1 = atom(in1)
      do j=1,3
        xca1(j)=xc(in1,j)
      enddo
      rdum = 2.0d0
121   continue       
      call FindNei( k, rdum )  
      if(inei.ge.3) then
        do i=3,inei
          ij = iatn(i)
          call CalcAng( xc(in1,1), xc(in1,2), xc(in1,3),
     *                  xc(k,1),   xc(k,2),   xc(k,3),
     *                  xc(ij,1),  xc(ij,2),  xc(ij,3),ang )
          if(ang.le.160.d0) then
            call PrtAng(nout,k,in1,ij,ang)     
            at2=atom(ij)
            do j=1,3
              xca2(j)=xc(ij,j)
            enddo
            goto 71
          endif          
        enddo          
      endif
c.....if no neighbors found within ~10 Angs - make dummy at 0,0,0               
      if(rdum.ge.10.0) then
        write(   *,384) atom(i1st), rdum
        write(nout,384) atom(i1st), rdum
        write(   *,386)
        write(nout,386)
        call AddDummyAtom(ifra,0.,0.,0.,idum,xca2(1),xca2(2),xca2(3))
        at2 = dlab(idum)
        goto 71
      endif
      rdum = rdum + 0.5D0 ! if no neigbours within rdum are found - increase it
      goto 121
c
c ================================
c  3 neighbors - sp2 C, N, O etc
c ================================
c      
30    continue
c.....decisions are based on whether the atomic group is planar or not...
      if( lPG ) then
c.......planar group
        if( is.eq.0 ) goto 300 
        if( is.eq.1 ) goto 301
        if( is.eq.3 ) goto 302
      else
c.......NON-planar group
        if( is.eq.0 ) goto 310 
        if( is.eq.1 ) goto 311
        if( is.eq.3 ) goto 312
      endif  
      line = '  '    
      write(line,399) is
399   format(' is = ',i4,' is not valid...')
      call Err('DefLoc',line(:ils(line)))
c
350   format(/' Label = ',i5)      
c
c     Planar: All 3 neighbors are different: m 
c     ----------------------------------------
c     m plane is in the plane of the atomic group
c
300   continue
      write(nout,350) 300
      ax1 = 'X'
      ax2 = 'Y' 
      asym(k) = 'm' 
      at1 = atom(isn(1))
      at2 = atom(isn(2))
      do j=1,3
        xca1(j) = xc(isn(1),j)
        xca2(j) = xc(isn(2),j)
      enddo
      goto 71
c
c     Planar : All 3 neighbors are the same: 3m
c     ------------------------------------------
c
302   continue
      write(nout,350) 303
      at1 = atom(i1st)
      at2 = atom(i2nd)
      ax1 = 'X'
      ax2 = 'Y' 
      asym(k) = '3m' 
      do j=1,3
        xca1(j) = xc(i1st,j)
        xca2(j) = xc(i2nd,j)
      enddo
      goto 71
c
c     Planar: Case 2 + 1: mm2 
c     -----------------------
c     '2' axis relates two 'same atoms' and goes thru 'different' atom
c
301   continue
      write(nout,350) 301
      ax1 = 'Z'
      ax2 = 'X' 
      asym(k) = 'mm2'
c     find atom which is 'different'
      do i=1,nnei(k)
        iat = jnei(k,i)
        if( iat.ne.isame(1,1) .and. iat.ne.isame(1,2) ) ia1 = iat ! atom which is different
      enddo         
c.....1 H-atom: Z-along midpoint of 2 other neighbors; X-along one of 'same' atoms
      if( itype(ia1).eq.1 ) then        
        i1 = isame(1,1)
        i2 = isame(1,2)
        call MidPoint( xc(i1,1),xc(i1,2),xc(i1,3),
     *                 xc(i2,1),xc(i2,2),xc(i2,3),
     *                    xm(1),   xm(2),   xm(3) ) 
        call AddDummyAtom( 0, xm(1),xm(2),xm(3), idum, 
     *                     xca1(1),xca1(2),xca1(3) )
        at1 = dlab(idum)          
        at2 = atom( i1 )
        do l=1,3
          xca2(l)=xc(i1,l)
        enddo
        goto 71
c.....2 H-atoms or no H-atoms: 
      else
        at1 = atom(ia1)
        at2 = atom(isame(1,1))
        do j=1,3
          xca1(j) = xc(ia1,j)
          xca2(j) = xc(isame(1,1),j)
        enddo
        goto 71
      endif
c
c     NON-Planar: All 3 neighbors are different: 1
c     --------------------------------------------
310   continue
      write(nout,350) 310
      ax1 = 'X'
      ax2 = 'Y' 
      asym(k) = 'NO' 
      at1 = atom(isn(1))
      at2 = atom(isn(2))
      do j=1,3
        xca1(j) = xc(isn(1),j)
        xca2(j) = xc(isn(2),j)
      enddo
      goto 71
c
c     NON-Planar: Case 2 + 1: m 
c     -------------------------
c     m relates the two 'same' atoms
c
311   continue
      write(nout,350) 311
      ax1 = 'X'
      ax2 = 'Y' 
      asym(k) = 'm'
      do i=1,nnei(k)
        iat = jnei(k,i)
        if( iat.ne.isame(1,1) .and. iat.ne.isame(1,2) ) ia1 = iat ! atom which is different
      enddo         
      i1 = isame(1,1)
      i2 = isame(1,2)
      call MidPoint( xc(i1,1),xc(i1,2),xc(i1,3),
     *               xc(i2,1),xc(i2,2),xc(i2,3),
     *                  xm(1),   xm(2),   xm(3) ) 
      call AddDummyAtom( 0, xm(1),xm(2),xm(3), idum, 
     *                   xca1(1),xca1(2),xca1(3) )
      at1=dlab(Idum)          
      at2=atom( ia1 )
      do j=1,3
        xca2(j)=xc(ia1,j)
      enddo
      goto 71
c
c     NON-Planar: All 3 neighbors are the same : 3m
c     ---------------------------------------------
c.....
312   continue
      write(nout,350) 312
      ax1 = 'Z'
      ax2 = 'X' 
      asym(k) = '3m'
c     get the center of the "ring" defined by 3 equivalent neighbors
      do l=1,3
        xm(l) = zero
        do i=1,3
          iat = jnei(k,i)
          if(ifra.eq.1) xm(l) = xm(l) + xf(iat,l)
          if(ifra.eq.0) xm(l) = xm(l) + xc(iat,l)
        enddo
        xm(l) = xm(l) / 3.00D00
      enddo        
      call AddDummyAtom( ifra, xm(1),xm(2),xm(3), idum, 
     *                   xca1(1),xca1(2),xca1(3) )
      at1 = dlab(idum)          
c     second axis - 1st nei (they are equivalent anyway)      
      at2 = atom(i1st)
      do l=1,3
        xca2(l) = xc(i1st,l)
      enddo
      goto 71
c
c ============================
c  4 neighbors - sp3 C, N etc
c ============================
c      
40    continue
c      write(*,*) ' 4 neighbors'
c      pause
c.....make decision based on the number of quivalent atom pairs 'is'
      goto(401,402,403,404,404,404) is
c.....is=0: all 4 atoms are UNequivalent - no symmetry
      ax1 = 'X'
      ax2 = 'Y'
      asym(k) = 'NO'
      ij = 0
      at1 = atom(isn(1))
      at2 = atom(isn(2))
      do j=1,3
        xca1(j) = xc(isn(1),j)
        xca2(j) = xc(isn(2),j)
      enddo
      goto 71
c.....is=1: 2 + 1 + 1 - symmetry m
401   continue
      ax1 = 'X'
      ax2 = 'Y'
      asym(k) = 'm'
      ij = 0
      id1 = isame(1,1)
      id2 = isame(1,2)
      do i=1,nnei(k)
        iat = isn(i)
        if( iat.ne.id1 .and. iat.ne.id2 ) then
          ij = ij+1
          if( ij.eq.1 ) i1 = iat        
          if( ij.eq.2 ) i2 = iat
        endif 
      enddo  
      at1 = atom(i1)
      at2 = atom(i2)
      do j=1,3
        xca1(j) = xc(i1,j)
        xca2(j) = xc(i2,j)
      enddo
      goto 71
c.....is=2: 2 + 2 - symmetry mm2 - OK?
402   continue 
      ax1 = 'Z'
      ax2 = 'Y'
      asym(k) = 'mm2'
      do i=1,is
        iat = isame(i,1)
        itypi = itype(iat)
        if(itypi.ne.1) then
          il1 = isame(i,1)
          il2 = isame(i,2)
          call MidPoint(xc(il1,1), xc(il1,2), xc(il1,3), 
     *                  xc(il2,1), xc(il2,2), xc(il2,3),
     *                      xm(1),     xm(2),     xm(3) )
          call AddDummyAtom( 0, xm(1),  xm(2),  xm(3), idum,
     *                       xca1(1),xca1(2),xca1(3) )
          at1 = dlab(idum)
          is1 = i
c          write(*,*) ' il1,il2,is1 = ',il1,il2,is1
c          write(*,'(3f10.4)') (xf(il1,l),l=1,3)
c          write(*,'(3f10.4)') (xf(il2,l),l=1,3)
c          write(*,'(3f10.4)') (xdum(ndum,l),l=1,3)
          goto 452
        endif
      enddo    
452   continue      
      do i=1,is
        if(i.ne.is1) then
          at2 = atom(isame(i,1))
          do j=1,3
            xca2(j)=xc( isame(i,1),j )
          enddo
          goto 71
        endif
      enddo
      call Err('DefLoc','4 neighbors: 2 + 2 case')
c.... is=3: 3 + 1 - symmetry 3m - OK
403   continue
      ax1 = 'Z'
      ax2 = 'X'
      asym(k) = '3m'
      at2 = atom( isame(1,1) )
      do j=1,3
        xca2(j) = xc( isame(1,1),j )
      enddo
c     choose the one which is NOT equiv
      do i=1,nnei(k)
        iat = jnei(k,i)
        im = 0
        do m=1,is
          if( isame(m,1).eq.iat .or. isame(m,2).eq.iat ) im = 1
        enddo      
        if(im.eq.0) then
          at1 = atom(iat)
          do j=1,3
            xca1(j)=xc(iat,j)
          enddo
        endif   
      enddo  
      goto 71
c.....is=6: all 4 atoms are equivalent - -42m (though, should be cubic) - OK?
404   continue
      ax1 = 'Z'
      ax2 = 'Y'
      asym(k) = '-42m'
      call MidPoint( xc(1,1),xc(1,2),xc(1,3),
     *               xc(2,1),xc(2,2),xc(2,3),
     *                 xm(1),  xm(2),  xm(3) ) 
      call AddDummyAtom( 0, xm(1),xm(2),xm(3), idum, 
     *                   xca1(1),xca1(2),xca1(3) )
      at1 = dlab(idum)
      at2 = atom(3)
      do j=1,3
        xca2(j) = xc(3,j)
      enddo
      goto 71
c
c ===========================================================
c  Almost done, but in order to get kappa sets properly we 
c  have to check all previous atoms 
c  ALSO DETERMINE EXCEPTIONS!!!!
c ===========================================================
c
71    continue
c
c     just want to create XD files w/o DB
c
      if( lsx ) then
        kset(k) = isft(k)
        itkset = max( itkset, kset(k) )
        chemcon = .false.
        goto 74
      endif
c
c     one of neighbors is disordered but the atom itself is not
c
      if( lndis .and. (.not.lkdis) ) then
        chemcon = .false.
c       let's find kappa set which has correct SCAT & 1's everywhere
        do 824 i=1,nkap
          if( nkapsf(i).ne.isft(k) ) goto 824
          if( abs( rkappa(i,1)-one ) .gt. 1.D-4 ) goto 824
          if( abs( rkappa(i,2)-one ) .gt. 1.D-4 ) goto 824
c         ok, looks like we found suitable kappa set
          write(nout,'(/a,i4)') ' Using old [IAM] kappa set = ',i
          kset(k) = i
          goto 74
824     continue
c       suitable kappa set was not found - let's add a new one
        goto 72
      endif
c
c     TODO : This is ugly, but I need to do it because if at least one
c     symmetry-generated atom was used in the defintion of the local 
c     coordinate system, then we are screwed... So we check coordinates
c     of each of atoms used in defition of loc.sys. against all coordinates
c     of symmetry-generated. If found a match, replace it with a dummy atom
c
      if( ia.ne.natoms ) call CheckSymEq(xca1,xca2,at1,at2,ax1,ax2)
c
c.....EXCEPTIONS
c
      write(nout,'(/a)') ' Calling DetExcept...'
      call DetExcept(k)  
      call upcase(group(k))  
      ij1 = ils(atom(k))
      ij2 = ils(group(k))
      if(iexp(k).eq.1) write(nout,4569) atom(k)(:ij1), group(k)(:ij2)
4569  format(/' Exception identified for ',a,
     *' : Functional group = [',a,']')      
c
c.....KAPPA SETS         
c
      if(k.eq.1) goto 72 ! first atom - first kappa set
      write(nout,'(/a)') ' Trying to determine the KAPPA SET : '
      do 710 i=1,k-1
c
c        H-atom
c       
         if( itype(k).eq.1 ) then
c          do NOT apply chemcon 
           if( .not. lhyd ) goto 72                    
c          tests for chemical equivalency
           if( isft(i).eq.isft(k)   .and.
     *          ihp(i).eq.ihp(k)    .and.   isp(i).eq.isp(k) .and. 
     *         ihnei(i).eq.ihnei(k) .and. group(i).eq.group(k) ) goto 73
c
c        "HEAVY" atoms
c
         else
c          do NOT apply chemcon 
           if( .not. lheav ) goto 72                    
c          equvalent atoms must have the same chem. environment & symmetry
c          we should NOT use 'iSimple' here - it fucks up kappa sets
           if( SameAtoms(k,i,1,1,0,1) ) goto 73
         endif      
c
710     enddo
c
c     new kappa set
c
72    itkset = itkset+1
      write(nout,'(/a,i4)') ' New kappa set = ',itkset
      kset(k) = itkset
      chemcon = .false.
      goto 74
c
c     old kappa set
c
73    kset(k) = kset(i)      
      write(nout,'(/a,i4)') ' Old kappa set = ',kset(i)
      chemcon = .true.
      ichem = i
c
C     OK, NOW WRITE DOWN COORDINATE SYSTEM DEFINITION
c
74    continue
cc
        atom1=atom(k)
        call erbras(atom1)
        atom2=atom(k)
        call ernum(atom2)
cc               
      ichemcon(k) = chemcon
      at = ' '
      at = atom(k)
      write(nout,'(/3a)') ' Symmetry = [',asym(k),']'
      line = ' '
      if(chemcon) then
        write(line,80) at,at1,ax1,at,at2,ax2,iut(k),isft(k),kset(k),
     *                 lmx(k),asym(k),atom(ichem) 
        ato1=at1
        call erUM(ato1)
        call erbras(ato1)
        ato2=at2
        call erUM(ato2)
        call erbras(ato2)
        ato3=atom(ichem)
        call erbras(ato3)
        ato=at
        call erbras(ato)
        ato4=at
        call ernum(ato4)
        write(dstring,'(I5)') kset(k)
        if(ato4.eq.'H') then
        write(line2,81) ax1,ax2,ato1,ato2,' QUA ',
     *  ' K'//adjustl(trim(dstring)), 'V0  M0  Q0'         
        else
        write(line2,81) ax1,ax2,ato1,ato2,' HEX ',
     *  ' K'//adjustl(trim(dstring)), 'V0  M0  Q0' 
        endif  
        if(resi11) then   
        write(line3,'(A,1X,A,2X,A)') '! CONPVM ',
     *  RESI2(k,1),RESI2(i,1)         
        else   
        write(line3,'(A,1X,A,2X,A,1X,A,1X,A)') '! CONPVM ',
     *  ato,'1',ato3,'1' 
        endif          
        jk = ils(line3)
        MCON(k,1)=line3(:jk)            
      else
        ato1=at1
        call erUM(ato1)
        call erbras(ato1)
        ato2=at2
        call erUM(ato2)
        call erbras(ato2)
        write(dstring,'(I5)') kset(k)
        write(line,80) at,at1,ax1,at,at2,ax2,iut(k),isft(k),kset(k),
     *                 lmx(k),asym(k)
        ato4=at
        call ernum(ato4)
        write(dstring,'(I5)') kset(k)
        if(ato4.eq.'H') then
        write(line2,81) ax1,ax2,ato1,ato2,' QUA ',
     *  ' K'//adjustl(trim(dstring)), 'V0  M0  Q0'         
        else
        write(line2,81) ax1,ax2,ato1,ato2,' HEX ',
     *  ' K'//adjustl(trim(dstring)), 'V0  M0  Q0' 
        endif        
      endif
80    format(a,1x,a,2x,a,2x,a,1x,a,1x,a,3x,'R ',2i3,2i4,3x,a,3x,a)  
81    format(A,A,3X,A,3X,A,3X,A,3X,A,3X,A)            
      ijk = ils(line)
      MOPRO(k,1)=line2(:ijk)
      write(   *,99) ' ', line(:ijk)
      write(imas,99) line(:ijk)
      write(nout,'(/a/)') ' Local coordinate system definition : '
      write(nout,99) ' ', line(:ijk)
cccccc
      if(asym(k).ne.'NO') then
        if(asym(k).eq. 'm') then
          if((ax1.eq.'X').and.(ax2.eq.'Y')) sym='mz'
          if((ax1.eq.'Y').and.(ax2.eq.'Z')) sym='mx'
          if((ax1.eq.'X').and.(ax2.eq.'Z')) sym='my'
          if((ax1.eq.'Y').and.(ax2.eq.'X')) sym='mz'
          if((ax1.eq.'Z').and.(ax2.eq.'Y')) sym='mx'
          if((ax1.eq.'Z').and.(ax2.eq.'X')) sym='my'
        endif
        if(asym(k).eq.'mm2') sym='mxmz'
        if(asym(k).eq.'3m') sym='3m'
        if(asym(k).eq.'cyl') sym='cy'
        ato=at
        call erbras(ato)
        if(resi11) then
      write(443,'(A,1X,A6,1X,A)') 'SYMPLM ',sym, RESI2(k,1)  
        else
      write(443,'(A,1X,A,1X,A,1X,A)') 'SYMPLM ',sym,ato,'1'
        endif
      endif   
!rkappa(nat,2),nkapsf(nat)      
C
C     VISUALIZATION OF LOCAL COORDINATE SYSTEM IN PLATON
C
      call Loc2Plat(k,xca1,xca2,ax1,ax2)
c
c     READING PSEUDOATOM DATABANK
c
      lreaddb = .false.
c     read db only in the following cases :
      if( ldb ) then
c       if k-th atom and its neighbors are ALL ordered
        if( .not.lkdis .and. .not.lndis ) lreaddb = .true.
c       if k-th atom is disordered & special request was made
        if( lkdis .and. ldd ) lreaddb = .true.
      endif      
      i9 = ils(atom(k))
      if(      lreaddb ) write(nout,414) atom(k)(1:i9)
      if( .not.lreaddb ) write(nout,415) atom(k)(1:i9)
414   format(/' Atom ',a,' :: will now read databank')
415   format(/' Atom ',a,' :: will skip reading databank')
c      
      if( lreaddb ) then
c
        call ReadDB(k)
c
        if(k.eq.1) then
          nkap=1
          nkapsf(nkap)=isft(k)
        endif   
        if( k.ne.1. and. kset(k).gt.nkap ) then
          nkap = kset(k)
          nkapsf(nkap) = isft(k)
        endif  
        rkappa(kset(k),1)=rkap1(1)
        rkappa(kset(k),2)=rkap1(2)
        do i=1,ia
          if(at1.eq.atom(i)) iar(k,3)=i
          if(at2.eq.atom(i)) iar(k,5)=i
        enddo
        Pv(k) = Pv1*occ(k)
        sPv(k) = sPv1
        do l=0,4
          do m=-l,l
            Plm(k,l,m) = Plm1(l,m)*occ(k)
          enddo
        enddo  
c
c     NO reading DB was requested or 
c     just want to create XD files w/o DB or
c     one of the neighbors is disordered but the atom itself is not
c
      else
        Pv(k) = Pval( itype(k) )*occ(k)   
        nkap = max(nkap,kset(k))
        rkappa(kset(k),1) = one
        rkappa(kset(k),2) = one
        nkapsf(kset(k)) = isft(k)
      endif
c
      write(nout,393) k, kset(k), nkapsf(kset(k))
393   format(/' CKS :: Atom=',i6,2x,'kset=',i4,2x,'nsft=',i4)
c
c     for both cases
c     
      do i=1,3
        if(ax1.eq.axes(i)) iar(k,1)=i
        if(ax2.eq.axes(i)) iar(k,2)=i
      enddo
      iar(k,4)=k
      do i=1,ia
        if(at1.eq.atom(i)) iar(k,3)=i
        if(at2.eq.atom(i)) iar(k,5)=i
      enddo
      do i=1,ndum
        if(at1.eq.dlab(i)) iar(k,3)=i+ia
        if(at2.eq.dlab(i)) iar(k,5)=i+ia 
      enddo
      iar(k,6)=iut(k)  ! jtf    - the order of displacement tensor
      iar(k,7)=isft(k) ! itbl   - scattering factor number
      iar(k,8)=kset(k) ! isfz   - kappa set
      iar(k,9)=4       ! lmax   - max. level of multipole expansion   
      iar(k,10)=1      ! isym   - site symmetry code
      iar(k,11)=0      ! ichcon - chemical constraint
      if(chemcon) iar(k,11)=ichem
c
99    format(100a)
      return
      end
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c prints environment of i-th atom
      subroutine PrtAtomNei(i,io)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      common /verbosity/ iverbose
c
      i1 = ils(atom(i))
c     neighbors
      write(io,500)  atom(i)(:i1), nnei(i)
500   format(/1x,'Covalently-bonded neighbors to ',a,' = ',i4,/)     
      write(io,502)
502   format(5x,'i',7x,'n',5x,'r(Ang)',4x,'atom',7x,'occ')   
      do j=1,nnei(i)
        jat = jnei(i,j)
        write(io,503) j, jat,  rnei(i,j), atom(jat), occ(jat)
      enddo
c
      if( iverbose.le.1 ) return
c     neighbor types     
      write(io,504) jtnei(i)
      do j=1,jtnei(i)
        write(io,*) ' '
        write(io,505) 'Type = ',j    
        write(io,505) '  atomic number = ',janei(i,j)
        jdum = jnnei(i,j)
        write(io,505) '  neighbors of this type = ',jdum
        write(io,506) ( atom(jsnei(i,j,k)),k=1,jdum)
      enddo
503   format(i6,i8,2x,f9.3,4x,a,2x,f7.4)
504   format(/4x,'There is(are) ',i4,' atom type(s) among neighbors')
505   format(4x,a,i4)
506   format(6x,'atoms : ', 6(a,1x), 100(/,14x,6(a,1x)) )
      return
      end 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for the case when choice of local symmetry is arbitrary, i.e. 
c mostly when local symmetry is '1' (or 'm' in case of planar 
c sp2-hybridization), copy sorted neighbors (see FindAllNeighbors)
c to array 'isn'  
c this is necessary because bloody structures are not always good
c [or because chemistry is crap!]
c
c so array 'isn' will contain neighbors sorted in descending 
c order by atomic number, then for the same atomic number 
c they are listed in ascending order in terms of hybridization,
c and finally for the same hybridization the atoms are listed
c in ascending order of the distance to k-th atom
c
      subroutine SortNei(k)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /isn/ isn(maxnei)
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      l = 0
      do i=jtnei(k),1,-1
        do j=1,jnnei(k,i)
          l = l + 1
          isn(l) = jsnei(k,i,j)
c          if( isn(l).eq.k ) then
c            write(*,*) ' k,i,j,jsnei(k,i,j) = ',k,i,j,jsnei(k,i,j)
c            stop 
c          endif
        enddo
      enddo
      if( l.ne.nnei(k) ) call Err('SortNei','l.ne.nnei(k)')
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c determines the number of neighbors of k-th atom and 
c sorts them on distance in ascending order
c also determines the number of atom types based on atomic numbers 
c 
c 081505 rewritten
c 093003 merged with subroutine NeiTyp 
c 011602 also determines among the neighbors the number of atom types
c        and number of neighbors of each type (type means atomic number)
c
c when determine neighbors there are 2 choices in terms of distances
c rmax <= 0 : determine covalently-bonded neighbors within 'tol1'
c       > 0 : determine neigbors within 'rmax+tol1' distance...
c
      subroutine FindNei(k,rmax)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
c.....neighbor info from FindNei
      dimension iatn(maxnei), rnd(maxnei) 
      dimension ndt(maxnei), natnum(maxnei), ndtnum(maxnei,maxnei)
      common /neig/ inei, iatn, rnd, ity, ndt, natnum, ndtnum
c.....
      logical la, ldis        
      dimension idum(maxnei)
      common /files/ imas, nout, ixp, idb, icon 
      common /iprint/ iprint  
      common /ftype/ iform
      common /verbosity/ iverbose
      data zero /0.00D00/
c
c.....generate neighbors & sort them based on distance
      call GetNei(k,rmax,ldis)
c.....print out info (if requested)
      if(iprint.eq.1) then      
        i1 = ils(atom(k))
        if( rmax.le.zero ) then
          write(nout,500) atom(k)(:i1), inei
        else
          write(nout,501) atom(k)(:i1), rmax, inei
        endif  
        write(nout,502)
        do j=1,inei
          jat = iatn(j)
          write(nout,503) j, iatn(j),  rnd(j), atom(jat), occ(jat)
        enddo
        if(ldis) write(nout,504) 
      endif
500   format(/1x,'Covalently-bonded neighbors to ',a,' = ',i4,/)     
501   format(/1x,'Neighbors of ',a,' within ',f7.3,' Ang = ',i4,/)     
502   format(5x,'i',7x,'n',5x,'r(Ang)',4x,'atom',7x,'occ')   
503   format(i6,i8,2x,f9.3,4x,a,2x,f7.4)
504   format(/4x,'Note : at least 1 neighbor is disordered')
c
c     figure out how many atoms types are among the nei and how 
c     many nei of each type, then sort the list of atom types by 
c     the atomic number
c
      ity = 0      
      do 511 i=1,inei
        iat = iatn(i)
        itypi = itype( iat )
c.......first neighbor - first type      
        if(i.eq.1) then
          ity=1
          ndt(ity)=1
          ndtnum(ity,1)=iat
          natnum(ity)=itypi
c.......other neighbors
        else
          do j=1,ity
c...........already exists
            if(itypi.eq.natnum(j)) then
              ndt(j)=ndt(j)+1      ! number of atoms of ity-th type
              ndtnum(j,ndt(j))=iat ! seq numbers of neighbors of ity-th type in /atoms/
              goto 510
            endif  
          enddo     
c.........new neigbour type
          ity=ity+1           
          ndt(ity)=1
          ndtnum(ity,1)=iatn(i)
          natnum(ity)=itypi
510       continue         
        endif
511   continue
c.....printout before sorting
c      if(iprint.eq.1) then
c        write(nout,'(/4x,a,i4,a)') 'There is(are) ',ity,
c     *               ' atom type(s) among neighbors BEFORE SORT'
c        do i=1,ity
c          write(nout,'(/4x,a,i4)')'Type = ',i    
c          write(nout,'( 4x,a,i4)')'  atomic number = ',natnum(i)
c          write(nout,'( 4x,a,i4)')'  neighbors of this type = ',ndt(i)
c          write(nout,'( 4x,a,10(a5,1x))')'  atoms : ', 
c     *                ( atnei(ndtnum(i,j)),j=1,ndt(i))
c        enddo
c      endif
c.....sort nei types by atomic number
520   la=.false.
      do 5 i=1,ity-1
        if( natnum(i+1).lt.natnum(i) ) then
c          exchange atom types
           ijunk=natnum(i)
           natnum(i)=natnum(i+1)
           natnum(i+1)=ijunk
c          exchange atom pointers
           do j=1,ndt(i)
             idum(j)=ndtnum(i,j) 
           enddo 
           do j=1,ndt(i+1)
             ndtnum(i,j)=ndtnum(i+1,j)
           enddo 
           do j=1,ndt(i)
             ndtnum(i+1,j)=idum(j)
           enddo 
c          exchange number of nei of this type
           ijunk=ndt(i)
           ndt(i)=ndt(i+1)
           ndt(i+1)=ijunk
c
           la=.true.
         endif  
5     continue      
      if(la) goto 520               
c.....
      if( iprint.eq.1 .and. iverbose.ge.2 ) then
        write(nout,'(/4x,a,i4,a)') 'There is(are) ',ity,
     *               ' atom type(s) among neighbors'
        do i=1,ity
          write(nout,'(/4x,a,i4)')'Type = ',i    
          write(nout,'( 4x,a,i4)')'  atomic number = ',natnum(i)
          write(nout,'( 4x,a,i4)')'  neighbors of this type = ',ndt(i)
          write(nout,530) ( atom(ndtnum(i,j)),j=1,ndt(i))
        enddo
      endif
530   format( 6x,'atoms : ', 6(a,1x), 100(/,14x,6(a,1x)) )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c finds neighbors of k-th atom and sorts them by distance
c rmax <= 0 : get neigbors from connectivity matrix
c       > 0 : determine neigbors within 'rmax+tol1' distance...
c ldis = .true. if k-th atom is connected to at least 1 disordered atom
      subroutine GetNei(k,rmax,ldis)
      include 'lsdb.inc'
      logical la      
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      common /files/ imas, nout, ixp, idb, icon 
      common /tols/ tol1, tol2, tol3
      common /covrad/ covrad(54)
      logical ldis
c     common /neig/ contains the neighbor info of k-th atom :
c     inei - number of neighbors 
c     iatn - sequence numbers of neig. in /atoms/
c     rnd  - distances in Angstroms (SORTED IN ACSENDING ORDER)
c     ity  - number of neighbor types (i=1..ity)
c     ndt(i) - number of neighbors of i-th type
c     natnum(i) - atomic number of i-th neighbor type
c     ndtnum(i,k), k=1...ndt(i) - seq. number of k-th nei of i-th type in /atoms/      
      dimension iatn(maxnei), rnd(maxnei) 
      dimension ndt(maxnei), natnum(maxnei), ndtnum(maxnei,maxnei)
      common /neig/ inei, iatn, rnd, ity, ndt, natnum, ndtnum
c
      data zero / 0.00D00 /, samexyz / 1.00D-02 /, tooclose / 0.5D0 /
c
      call izero(iatn,maxnei)
      call rzero(rnd,maxnei)
      call izero(ndt,maxnei)
      call izero(natnum,maxnei)
      call izero(ndtnum,maxnei*maxnei)
c
      xk = xc(k,1)
      yk = xc(k,2)
      zk = xc(k,3)
      inei = 0
      itypk = itype(k)
      ldis = .false.
c
c     find neighbors
c
      do 10 i=1,natoms
        if(i.eq.k) goto 10
c       correct treatment of disorder of k-th atom (should be connected
c       either to ordered atom or disordered atom with the same occupation       
        if( iaf(k).eq.-1 .and. iaf(i).eq.-1 ) then
c         skip disordered atoms with different occupation
          docc = occ(k) - occ(i)
          if( abs(docc).gt.1.D-2 ) goto 10
        endif        
c       calculate distance
        r1 = Dist( xk, yk, zk, xc(i,1), xc(i,2), xc(i,3) )
c       check the desired cutoff   
        if( rmax.gt.zero  ) then
          if( r1.gt.rmax ) goto 10
        else
          rcov = covrad( itypk ) + covrad( itype(i) ) + tol1       
          if( r1.gt.rcov ) goto 10
        endif
c       make sure that the atom is not too close (for example when
c       disorder is present)
        if( r1.lt.tooclose ) goto 10
c       make sure that this neighbor does not occupy the same position
c       as one of the already accepted neighbors (disorder)
        do j=1,inei
          jn = iatn(j)
          rj = Dist( xc( i,1), xc( i,2), xc( i,3), 
     *               xc(jn,1), xc(jn,2), xc(jn,3) ) 
          if( rj.lt. samexyz ) goto 10 
        enddo         
c       Ok, fill in part of /neig/ arrays
        inei = inei + 1 
        if( inei.gt.maxnei ) then
           write(   *,999) k, atom(k), rmax, maxnei
           write(nout,999) k, atom(k), rmax, maxnei
           stop
        endif
        iatn(inei) = i      ! sequence number of inei-th neighbor in /atoms/
        rnd(inei) = r1      ! distance to inei-th neighbor
        if( iaf(i).eq.-1 ) ldis = .true.
10    continue
999   format(/' GetNei :: atom = ',i8,2x,a,/,
     *        '           rmax = ',f12.3,/,
     *        '           inei exceeds dimension maxnei ( ',i6,' )'/)                 
c
c     sort neighbors on distance 
c
12    la=.false.
      do 5 i=1,inei-1
        if(rnd(i+1).lt.rnd(i)) then
          call ExR(rnd(i),rnd(i+1))
          call ExI(iatn(i),iatn(i+1))
          la=.true.
        endif
5     continue
      if(la) goto 12
c
      return
      end      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ExR(a,b)
      c=a
      a=b
      b=c
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine ExI(i,j)
      k=i
      i=j
      j=k
      return
      end 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c called whenever a new dummy atom needs to be created 
c it also checks input coordinates of the "to-be-created" dummy
c againts all dummy atoms already in the list. 
c on the input {xf,yf,zf} are the input coordinates of the dummy 
c                         atom "to-be-created" (see 'ifl')
c              'ifl' - flag for the type of {xf,yf,zf} supplied
c                    = 1 : fractional
c                    = 0 : cartesian     
c on the output 'idum' is the serial number of the new or old dummy atom
c              {xc,yc,zc} are ALWAYS CARTESIAN coordinates of the 
c                         idum-th (requested) dummy atom
      subroutine AddDummyAtom(ifl,xi,yi,zi,idum,xc,yc,zc)
      include 'lsdb.inc'
      common /ifra/ ifra
      character dlab*6
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      data zero /0.00D00/
c
c     if we can make use of fractional coordinates
c
      if( ifra.eq.1 ) then
c       input in cartesian
        if( ifl.eq.0 ) then
c         save cartesian
          xc = xi
          yc = yi
          zc = zi 
c         get fractional
          call rabs(1,xc,yc,zc,xi,yi,zi)  
        endif
c.....only cartesian
        xc = xi
        yc = yi
        zc = zi 
      endif      
c
c.....first dummy      
      if( ndum.eq.0 ) then
        ndum = 1
        call DummyLab(ndum)
        xdum(ndum,1) = xi
        xdum(ndum,2) = yi
        xdum(ndum,3) = zi
c.....check ccordinates of dummies already in the list 
      else
c       loop over all dummies in the list
        do i=1,ndum
c         get cartesian coordinates in {xd,yd,zd}
          if( ifra.eq.1 ) then
            call rabs(0,xdum(i,1),xdum(i,2),xdum(i,3),xd,yd,zd)
          else
            xd = xdum(i,1) 
            yd = xdum(i,1) 
            zd = xdum(i,1) 
          endif  
c         calculate the distance between two dummies
          r = Dist( xc,yc,zc, xd,yd,zd )
          if( r.lt.1.D-01 ) then
            idum = i
c           return cartesian coordinates of idum-th dummy atom
            xc = xd
            yc = yd
            zc = zd
            return
          endif
        enddo
c       ok, no matching dummy found... make new dummy        
        ndum = ndum + 1
        call DummyLab(ndum)
        xdum(ndum,1) = xi
        xdum(ndum,2) = yi
        xdum(ndum,3) = zi
      endif  
      idum = ndum
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c label for i-th dummy atom - writes directly to common block /dummy/
      subroutine DummyLab(i)
      include 'lsdb.inc'
      character dlab*6  
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
c
      dlab(i)=' '
      if(i.le.9) then
        write(dlab(i),'(a3,i1)') 'DUM',i
      elseif(i.ge.10.and.i.le.99) then
        write(dlab(i),'(a3,i2)') 'DUM',i
      elseif(i.ge.100.and.i.le.999) then
        write(dlab(i),'(a3,i3)') 'DUM',i
      else
        call Err('DummyLab','DUMMY label exceeds format...')
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c puts all references to neighbors of k-atom into array 'id'
c the last atom in array 'id' is the k-th atom itself
c nn - final number of atoms in array 'id'
      subroutine CollectNei(k,id,nn)
      include 'lsdb.inc'
      dimension id(maxmem)
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each neighbor type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      call izero(id,maxmem)
      nn = 0
      do i=1,nnei(k)
        nn = nn + 1
        iat = jnei(k,i)
        id(nn) = iat
      enddo
      nn = nn + 1
      id(nn) = k
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c calculate angle between atoms 1-2-3 with cartesian 
c coordinates x1-x2-x3. returned angle is in degrees
      subroutine CalcAng(r1,r2,r3,r4,r5,r6,r7,r8,r9,ang)
      parameter( pi=3.1415926535897932384626433832795028841971694D0 )
      dimension x1(3), x2(3), x3(3)
      dimension d1(3), d2(3)
      data zero /0.00D00/
      ang=zero
      x1(1)=r1
      x1(2)=r2
      x1(3)=r3
      x2(1)=r4
      x2(2)=r5
      x2(3)=r6
      x3(1)=r7
      x3(2)=r8
      x3(3)=r9
      ang=zero
      dum1=zero
      dum2=zero
      dum3=zero
      do j=1,3
       d1(j) = x1(j) - x2(j)
       d2(j) = x3(j) - x2(j)   
       dum1 = dum1 + d1(j)**2 
       dum2 = dum2 + d2(j)**2 
       dum3 = dum3 + d1(j)*d2(j)            
      enddo
      dum1=sqrt(dum1)
      dum2=sqrt(dum2) 	   	
      ang=(180.d0/pi)*acos( dum3/(dum1*dum2))
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prints angle between three atoms to unit 'iout'
      subroutine PrtAng(io,i,j,k,ang)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      i1 = ils( atom(i) )
      j1 = ils( atom(j) )
      k1 = ils( atom(k) )
      write(io,837) atom(i)(:i1), atom(j)(:j1), atom(k)(:k1), ang
837   format(/' Angle   ',a,' - ',a,' - ',a,' = ',f8.2,' deg.')        
      return
      end      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c extends the X--H distance to neutron value where H atom is (k) and
c its parent X-atom is (l)
c actual group we have to consider = Y--X--H
c
      subroutine ExtendH(l,k)
      include 'lsdb.inc'
      dimension id(maxmem)
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
c     ihp - atomic number of parent atom 
c     isp - hybridization of parent atom (total number of nei)
c     ihnei - number of attached H-atoms to parent atom
      character*17, dimension(nat,1):: RESI1
      character*10, dimension(nat,1):: RESI2
      logical resi11
      common /resi1/ RESI1, RESI2, resi11
      common /forH/ isp(nat), ihp(nat), ihnei(nat)
      common /files/ imas, nout, ixp, idb, icon
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /ifra/ ifra   
      common /iprint/ iprint  
      logical led, ldd
      common /disorder/ led, ldd
      common /xhbonds/ irb, icb, ires
      character group*50, line*255, adum1*20, adum2*20, ax(3)*1
      character line2*100
      character line1*255
      character q1*10, q2*10
      character*100, dimension(nat,1)::RES
      common /cons/ RES
      real t
      logical laC, lc, lp, lPG
      data zero /0.00D00/, ii /0/, line /' '/
      data ax / 'X','Y','Z' /
      save ii, line
c
      iprint = 0
c      call FindNei( l, zero )
      rnew = zero     
      group = ' '  
      l1 = ils(atom(l))
      k1 = ils(atom(k))
c     if H-atom is disordered and request was not to extend H-distance
      if( iaf(k).eq.-1 .and. .not.led ) then
        write(nout,386)atom(l)(:l1),atom(k)(:k1),'H-ATOM IS DISORDERED'
        return
      endif
c
      write(nout,'(/3a)') ' Trying to extend hydrogen ',
     *                    atom(k)(:k1),' to standard distance :'
c     parent (X) atom info
      jp = ihp(k)        ! parent atom atomic number
      jsp = isp(k)       ! hybrydization of parent atom  
      jhnei = ihnei(k)   ! number of H-atom connected to parent atom  
c     make choice base on parent (X)-atom
      if( jp.eq. 6 ) goto  60
      if( jp.eq. 7 ) goto  70
      if( jp.eq. 8 ) goto  80     
      if( jp.eq.16 ) goto 160     
      write(nout,386)atom(l)(:l1),atom(k)(:k1),
     *  'DON''T KNOW HOW TO HANDLE'
      return
c
c____ C -- H ___________________________________________________________
cKNJmodyf
60    continue
c     first of all, let's decide if all atoms connected to parent (X)
c     C-atom are carbons or NOT ,i.e. what are the 'Y'-atoms        
      laC = .true.
c     loop over types  
      do i=1,nnei(l)
        iat = jnei(l,i)
        itypi = itype(iat) 
        if( itypi.ne.6 .and. itypi.ne.1 ) then
          laC = .false.
          goto 61
        endif
      enddo
61    continue
c
c     sp3 carbon
c
      if( jsp.eq. 4 ) then
c        methyl 
         if( jhnei.eq.3 ) then
            if(laC) then
              rnew = 1.077D00
              group = 'C-C-H3' 
            else
              rnew = 1.077D00
              group = 'X-C-H3' 
            endif
         endif 
c        primary
         if( jhnei.eq.2 ) then
            if(laC) then
              rnew = 1.092D00
              group = 'C2-C-H2' 
            else
              rnew = 1.091D00
              group = 'X2-C-H2' 
            endif
         endif 
c        secondary
         if( jhnei.eq.1 ) then
            if(laC) then
            rnew = 1.099D00
              group = 'C3-C-H' 
            else
            rnew = 1.098D00
              group = 'X3-C-H' 
            endif
         endif 
      endif
c
c     sp2 carbon
c
      if( jsp.eq. 3 ) then
        call IsItCyclic(l,lc,lp,lx,jr,mem)
c       aromatic
        if( jhnei.eq.1 .and. lc .and. lp ) then      
          rnew = 1.083D00
          group = 'aromatic C-H'
        else
          rnew = 1.082D00
          group = 'C-C=C-H'
        endif
      endif      
      goto 110
c
c
c     sp carbon
c
      if( jsp.eq. 2 ) then
          rnew = 1.042D00
          group = 'C=Csp-H'
      endif      
      goto 110
c
cKNJstart
c____ N -- H ___________________________________________________________
c
70    continue
c     X2-N-H
      if( jsp.eq.3 .and. jhnei.eq.1 ) goto 71
c     X-N-H2 planar or not?  
      if( jsp.eq.3 .and. jhnei.eq.2 ) then
        call CollectNei(l,id,imem)
        call Planarity(0,lPG,imem,id,3)
        if( lPG ) goto 72
        if(.not. lPG) goto 73
      endif
c
c     X-N-H3
      if( jsp.eq.4 ) then
        rnew = 1.036D00
        group = 'N+-H'
      endif
c
      goto 110
c
71     continue 
c     X2-N-H
        iC = 0
        do i=1,nnei(l)
          iat = jnei(l,i)
          itypi = itype(iat)
          in = nnei(iat)
          if( itypi.eq.6 .and. in.eq.3 ) iC = iC + 1
        enddo
        if( iC.eq.2) then
          rnew = 1.030D00
          group = 'C2-N-H'
          goto 110
        else
          rnew = 1.027D00
          group = 'X2-N-H'
        endif
        goto 110
c
72      continue
c     X-N-H2 planar
        iC = 0
        do i=1,nnei(l)
          iat = jnei(l,i)
          itypi = itype(iat) 
          if( itypi.eq.6) iC = iat
        enddo
        if( iC.ne.0 .and. nnei(iC).eq.3 ) then
         call IsItCyclic(iC,lc,lp,lp2,ir,mem) 
         if( lc .and. lp ) then
          rnew = 1.010D00
          group = 'C(ar)-N-H2pl'
          goto 110
         else
          ict = 0
          do j=1,nnei(iC)
              jat = jnei(iC,j)
              itypj = itype(jat)
          if( itypj.eq.8 .and. nnei(jat).eq.1) ict = 1
          enddo
          endif
          if( ict.eq.1) then
          rnew = 1.010D00
          group = 'C-N-H2 amido'
          goto 110
          else 
          rnew = 1.012D00
          group = 'C-N-H2pl'
          endif
        if( iC.ne.0 .and. nnei(iC).eq.4 ) then
          rnew = 1.002D00
          group = 'Csp3-N-H2pl'
        endif
        if( iC.eq.0 ) then        
          rnew = 1.015D00
          group = 'X-N-H2pl'
        endif
        endif
        goto 110
c
73     continue
c      X-N-H2 not planar 
        iC = 0
        do i=1,nnei(l)
          iat = jnei(l,i)
          itypi = itype(iat) 
          if( itypi.eq.6) iC = iat
         enddo
        if( iC.ne.0 .and. nnei(iC).eq.3 ) then
         call IsItCyclic(iC,lc,lp,lp2,ir,mem) 
         if( lc .and. lp ) then
          rnew = 1.024D00
          group = 'C(ar)-N-H2py'
          goto 110
         else
          rnew = 1.015D00
          group = 'X-N-H2'
         endif
        endif
c
      goto 110
cKNJ
c____ O -- H ___________________________________________________________
cKNJ
80    continue
c
c     H2O :
c
      if( jsp.eq.2 .and. jhnei.eq.2 ) then
        rnew = 0.9584D00
        group = 'water'
        goto 110
      endif
c
c     Y-O-H :       
c
      if( jsp.eq.2 .and. jhnei.eq.1 ) then
c       get info about Y-atom, i.e. Y-O-H     
        iY = jnei(l,2)
c       find where the O-atom-type info is stored
        iot = 0
        do i=1,jtnei(iY)
          if( janei(iY,i).eq.8 ) iot = i
        enddo        
        if( iot.eq.0 ) call Err('ExtendH','Error 1 in Y-O-H')
c       only 1 oxygen around Y-atom        
       if (itype(iY) .eq. 6) goto 81
       if (itype(iY) .ne. 6) goto 82 
c
81     continue
c      Y = Carbon 
        if( jnnei(iY,iot).eq.1 ) then
         call IsItCyclic(iY,lc,lp,lx,jr,mem)
c        aromatic
         if( lc .and. lp ) then
          rnew = 0.992D00
          group = 'O-H in phenol'
          goto 110  
         endif
         if( nnei(iY).eq.4 ) then
          rnew = 0.970D00
          group = 'Csp3-O-H in alcohol'
          goto 110
         else 
          rnew = 0.980D00
          group = 'C-O-H any'
         endif 
        goto 110
        endif
c
c       2 oxygens around Y-atom
        if( jnnei(iY,iot).eq.2 ) then
c         find "other" O-atom     
          io2 = 0
          do i=1,nnei(iY)
            iat = jnei(iY,i)
            if( itype(iat).eq.8 .and. iat.ne.l ) io2 = iat
          enddo 
          if( io2.eq.0 ) call Err('ExtendH','Error 2 in Y-O-H')
c         now check environment of "other" oxygen
          if( nnei(io2).eq.1 ) then
            rnew = 1.018D00
            group = 'O-H in acid'
            goto 110
          else          
            rnew = 0.980D00
            group = 'O-H in ether'
          endif
        endif
       goto 110
cKNJ
cPMD start
c       3 or 4 oxygens around Y-atom 
c       (dirty way of defining p-O-H in phosphonic acid, phosphate etc.)      
82      continue
c       Y different than Carbon 
        if( jnnei(iY,iot).eq.3 .or. jnnei(iY,iot).eq.4 ) then
          rnew = 1.015D00
          group = 'O-H in Pacid'
          goto 110
        else
          rnew = 0.983D00
          group = 'X-OH'
        endif    
CPMD end
      endif 
      goto 110
cKNJ
c____ S -- H ___________________________________________________________
c
160   continue
      rnew = 1.338D00
      group = 'S-H'
      goto 110
cKNJ
c____ Ok, we've done all we can. Now let's sort it out _________________
c
110   continue
c
c     could not identify functional group...
c
      if( group.eq.' ' ) then
c        write(nout,'(3x,a)') 'Could not identify functional group :('      
        write(nout,386) atom(l)(:l1), atom(k)(:k1),'NOT IDENTIFIED'
c       but still want to use RESET BOND and CON instructions
c        return
c       calculate current distance and keep it 
        rnew = Dist( xc(l,1), xc(l,2), xc(l,3),  
     *               xc(k,1), xc(k,2), xc(k,3)  )
        goto 333
      endif   
c
c     actually extend the distance 
c
c      write(nout,'(3x,2a)') 'Functional group = ',group
      call ExtendBond( xc(l,1), xc(l,2), xc(l,3), 
     *                 xc(k,1), xc(k,2), xc(k,3), rnew, rold)      
      write(nout,386) atom(l)(:l1), atom(k)(:k1),  
     *                group(:ils(group)), ' =',rold, rnew
386   format(3x,'X-H group : ',a,'-',a,' >',2x,a,a,f7.3,' /',f7.3,
     * ' Angstroms')
c      write(nout,387) atom(l)(:l1), atom(k)(:k1), rold
c      write(nout,388) atom(l)(:l1), atom(k)(:k1), rnew
c387   format(3x,'Current  distance ',a,' - ',a,' = ',f8.3,' Angstroms') 
c388   format(3x,'Standard distance ',a,' - ',a,' = ',f8.3,' Angstroms') 
c
c     modify also fractional (if necessary), i.e. ifra=1
c
      if(ifra.eq.1)
     * call rabs(1,xc(k,1),xc(k,2),xc(k,3),xf(k,1),xf(k,2),xf(k,3))
c
c     repeat neighbor search for parent atom
c
      call FindNei( l, zero )
c
333   continue
c
c     write RESET BOND instruction to 'reset.bond' file
c
c     initialize RESET BOND line
c      write(*,*) ' ii = ',ii
      if( ii.eq.0 ) then
        line = 'RESET BOND  '
        line2 = 'URATIO  '
        ii = ils(line)+3
      endif  
c      write(*,'(3a)') '[',line,']'
      write(line(ii:),339) atom(l)(:l1), atom(k)(:k1), rnew
      t=0
      if(itype(l).eq.6)  t=1.2
      if(itype(l).eq.7)  t=1.5
      if(itype(l).eq.8)  t=1.5 
      if(itype(l).eq.16)  t=1.5
       q1=atom(l)
       call erbras(q1)
       q2=atom(k)
       call erbras(q2)
       if(resi11) then
      write(445,'(A,A,2X,A,1X,f6.3,2X,A)') 'DISTAN   ',  
     * RESI2(l,1),RESI2(k,1), rnew, '0.001'  
      write(line2(ii:),'(A,2X,A,1X,f6.1,2X,A)')  RESI2(l,1),
     * RESI2(k,1),t,'0.01'         
       else
      write(445,'(A,A5,1X,I1,2X,A5,1X,I1,1X,f6.3,1X,A)') 'DISTAN   ',  
     * q1(:l1),1,q2(:k1),1, rnew, '0.001'  
      write(line2(ii:),'(A5,1X,I1,2X,A5,1X,I1,1X,f6.1,1X,A)')  q1(:l1),
     * 1,q2(:k1),1,t,'0.01'  
       endif            
c      ii = ils(line) + 3
c      if( ii.ge.72 )  then
        write(irb,'(a)') line(:ils(line))
!        write(ires,'(a)') line2(:ils(line2))
        RES(k,1)= line2(:ils(line2))
!       write(445,'(A)') RES(k,1)
        ii = 0
c      endif  
339   format(a,' ',a,f6.3)
c
c     write CONSTRAIN bond instruction to 'const.bond' file
c
      write(icb,320) atom(l)(:l1), atom(k)(:k1), rnew
320   format('!',/,'! ',a,' - ',a,' = ',f9.3,' Angstroms',/,'!')      
      do i=1,3
        adum1 = ' '
        adum2 = ' '
        write(adum1,321) ax(i), l
        write(adum2,321) ax(i), k
        call RmBlanks( adum1 )
        call RmBlanks( adum2 )
        line1 = ' '
        write(line1,322) adum1(:ils(adum1)), adum2(:ils(adum2))
        write(icb,'(a)') line1(:ils(line))
      enddo  
321   format(a,'/',i8)
322   format('CON 1 ',a,' -1 ',a,' = 0')
c      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c extends the distance of atom with cartesian coodinates x2(1..3) along
c the bond x1(1..3)--->x2(1..3) at 'rnew' (in Angstroms)
      subroutine ExtendBond(x1,y1,z1,x2,y2,z2,rnew,rold)
      data zero / 0.00D00 /
c     calculate initial distance in Angstroms
      rold = zero
      dx = x2 - x1
      dy = y2 - y1
      dz = z2 - z1  
      rold = sqrt( dx*dx + dy*dy + dz*dz )
c     calculate delta
      delta = rnew / rold
c     extend cartesian
      x2 = delta * dx + x1        
      y2 = delta * dy + y1        
      z2 = delta * dz + z1        
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determine the coordinates of the midpoint between
c two atoms {x1,y1,z1} and {x2,y2,z2}. 
      subroutine MidPoint(x1,y1,z1,x2,y2,z2,xm,ym,zm)
      include 'lsdb.inc'
      xm = ( x1 + x2 ) / 2.D0
      ym = ( y1 + y2 ) / 2.D0
      zm = ( z1 + z2 ) / 2.D0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c checks if any of the atoms used in definition of local coordinate
c system is symmetry-generated.. if so, replace it(them) with dummy 
c atom(s)....
      subroutine CheckSymEq(xca1,xca2,at1,at2,ax1,ax2)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      character ax1, ax2, at1*8, at2*8
      character dlab*6
      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      common /files/ imas, nout, ixp, idb, icon
      common /cryst/ natoms, isymeq(nat), nsymeq(nat)
      dimension xca1(3), xca2(3)
      data thresh /1.00D-03/
      do i=ia+1,natoms
c       check 1st atom  
        r1 = Dist(xca1(1),xca1(2),xca1(3),xc(i,1),xc(i,2),xc(i,3))
        if( r1.le. thresh ) then
          ndum = ndum + 1
          call DummyLab(ndum)
          write(nout,100) at1(:ils(at1)), ax1, ndum
          at1 = dlab(ndum)
          do l=1,3
c            xdum(ndum,l) = xca1(l)
            xdum(ndum,l) = xf(i,l)
          enddo 
        endif 
c       check 2nd atom  
        r2 = Dist(xca2(1),xca2(2),xca2(3),xc(i,1),xc(i,2),xc(i,3))
        if( r2.le. thresh ) then
          ndum = ndum + 1
          call DummyLab(ndum)
          write(nout,100) at2(:ils(at2)), ax2, ndum
          at2 = dlab(ndum)
          do l=1,3
c            xdum(ndum,l) = xca2(l)
            xdum(ndum,l) = xf(i,l)
          enddo 
        endif 
      enddo
100   format(/' Atom ',a,' used in definition of ',a,'-axis is',
     *        ' symmetry-generated; replace with dummy ',i4)      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c determines exceptions (if applicable)
c
c***********************************************************************
      subroutine DetExcept(k)
      include 'lsdb.inc'
      character group(nat)*40
      common /except/ iexp(nat), group
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      common /files/ imas, nout, ixp, idb, icon 
      common /iprint/ iprint  
      common /isn/ isn(maxnei)
      logical lPG, lcyc, lPG2
      dimension id(maxmem)
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each neighbor type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      common /verbosity/ iverbose
c
c      write(*,*) ' k = ',k
96    format(/100a)
97    format(/' DetExcept :: checking ',a,' group...')
      iexp(k) = 0
      group(k) = ' '
      itypK = itype(k)
      if( itypK.eq.1 ) goto 1
      if( itypK.eq.6 ) goto 6
      if( itypK.eq.7 ) goto 7
      if( itypK.eq.8 ) goto 8
      if( itypK.eq.9 ) goto 9
      if( itypK.eq.15 ) goto 15
      if( itypK.eq.16 ) goto 16
      if( itypK.eq.17 ) goto 17
      if( itypK.eq.35 ) goto 35
      return
c
c---- HYDROGEN atom ------------------------------------------------------
c
1     continue
c     parent atom
      iP = jnei(k,1)
c     this is case when we have UNBONDED H-ATOM (can happen in fucking
c     macromolecular structures - they are so accurate... ) 
      if( iP.eq.0 ) return
c
      if( itype(iP).eq.8 .and. nnei(iP).eq.2 ) goto 1000
      if( itype(iP).eq.7 .and. nnei(iP).eq.3 ) goto 1010
      return
c
c     x-o-H
c      
1000  continue
      iexp(k) = 1
      write(nout,97) 'x-o-H'
c     check parent atom  
      if( jtnei(iP).le.1 .or. jtnei(iP).gt.2 ) then
        write(nout,999)  'x-o-H [1]'
        goto 111
      endif  
      iX = jsnei(iP,2,1)
cPMD start
c     x - is carbon 
      if( itype(iX).eq.6 ) goto 1001
c     x - is phosphorus
      if( itype(iX).eq.5 ) goto 1002
c
      if( itype(iX).eq.15 ) goto 1005
c     x - other than carbon or phosphorus
      if( itype(iX).ne.6 .or. itype(iX).eq.15 ) then   
        write(nout,999)  'x-o-H [2]'
        goto 111
      endif    
c     c-o-H in alcohols phenols or carboxylic acids     
1002  continue
      group(k) = 'b-o-H'
      return
      goto 111 
1001  continue
cPMD end       
      if( nnei(iX).eq.4 ) then
        group(k) = 'alcohol'
        return
      endif
      if( nnei(iX).eq.3 ) then
        call IsItCyclic(iX,lcyc,lPG,lPG2,iXr,mem)
        if( lcyc .and. lPG ) then
          group(k) = 'phenol'
          return
        endif  
        if( .not.lcyc .and. janei(iX,jtnei(iX)).eq.8 .and. 
     *                      jnnei(iX,jtnei(iX)).eq.2 ) then
          group(k) = 'carboxylic'
          return
        endif 
      endif
      write(nout,999)  'c-o-H'
      goto 111
      return   
cPMD  start
c     p-o-H in all 4-coordinated phosporus derivatives   
1005  continue
      if( nnei(iX).eq.4 ) then
          group(k) = 'phosp... acid'
          return
      endif
      write(nout,999)  'p-o-H'
      goto 111
      return        
cPMD
cPMD     x2-n-H
cPMD      
1010  continue
      iexp(k) = 1
      write(nout,97) 'x2-n-H'
      iN = jsnei(k,jtnei(k),1)
cPMD       check planarity of N-atom-group
        call CollectNei(iN,id,imem)
        idp = 0
        if( iverbose.ge.3 ) idp = 1
        call Planarity(idp,lPG,imem,id,3)
        if( lPG ) then
          call IsItCyclic(iN,lcyc,lPG, lPG2,iCr,mem) 
          if( lcyc .and. lPG ) then
          group(k) = 'aromatic'
          else
          group(k) = 'planar amine'
	  endif
        else
          group(k) = 'pyramidal amine'
        endif
        return        
      write(nout,999)  'x2-n-H'
      goto 111
      return        
cPMD end
c
c---- CARBON atom ------------------------------------------------------
c
6     continue
c 
      call IsItCyclic(k,lcyc,lPG,lPG2,iCr,mem)
cPMD  c2-C=o or c=Cc-o-x
      if( nnei(k).eq.3 .and. .not. ( lcyc .and. lPG ) .and.
     * janei(k,1).eq.6 .and. jnnei(k,1).eq.2 .and.
     * janei(k,2).eq.8 .and. jnnei(k,2).eq.1 ) goto 6050
c     c-C-o2
      if( nnei(k).eq.3 .and. 
     * janei(k,1).eq.6 .and. jnnei(k,1).eq.1 .and.
     * janei(k,2).eq.8 .and. jnnei(k,2).eq.2 ) goto 6000
cPMD  n-C-o2
      if( nnei(k).eq.3 .and. 
     * janei(k,1).eq.7 .and. jnnei(k,1).eq.1 .and.
     * janei(k,2).eq.8 .and. jnnei(k,2).eq.2 ) goto 6005
c     possible n-C-(c/h)3
      if( nnei(k).eq.4 .and. janei(k,jtnei(k)).eq.7 .and.
     *                       jnnei(k,jtnei(k)).eq.1 ) goto 6010 
c     C-o-X (C in ring, o-atom is not in ring)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,1).eq.6 .and. jnnei(k,1).eq.2 .and. 
     *   janei(k,2).eq.8 .and. jnnei(k,2).eq.1 ) goto 6020 
cKNJ  c2-C-b
      if( nnei(k).eq.3 .and. 
     * janei(k,1).eq.5 .and. jnnei(k,1).eq.1 .and.
     * janei(k,2).eq.6 .and. jnnei(k,2).eq.2 ) goto 6095     
cKNJ  C-n= (C in ring, 3 n as neighbours)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,jtnei(k)).eq.7 .and. jnnei(k,jtnei(k)).eq.3 ) goto 6060 
cPMD    C-n= (C in ring, 2 n as neighbours)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,jtnei(k)).eq.7 .and. jnnei(k,jtnei(k)).eq.2 ) goto 6035
cKNJ  n2-C-o (C in ring)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,1).eq.7 .and. jnnei(k,1).eq.2 .and. 
     *   janei(k,2).eq.8 .and. jnnei(k,2).eq.1 ) goto 6080       
cKNJ    n-C(=o)-c, n-C(-oh)-c, n-Cc-o (C in ring, o and n as neighbours)
      if( nnei(k).eq.3 .and. lcyc .and. lPG .and.
     *    janei(k,1).eq.6 .and. janei(k,2).eq.7 .and. janei(k,3).eq.8)
     *    goto 6070
cPMD    C-n= (C in ring, one n as a neighbor)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,jtnei(k)).eq.7 .and. jnnei(k,jtnei(k)).eq.1 ) goto 6030 
      if( nnei(k).eq.3 .and. lcyc .and. lPG2 .and.       
     *   janei(k,jtnei(k)).eq.7 .and. jnnei(k,jtnei(k)).eq.1 ) goto 6031 
c     c1=Cx1-c3 (where C, x1 & c3 are in the planar ring)     
      if( nnei(k).eq.3 .and. lcyc .and. lPG .and. 
     *   janei(k,jtnei(k)).eq.6 .and. jnnei(k,jtnei(k)).eq.2 ) goto 6040
cKNJ   C4-s-(?)
      if( nnei(k).eq.4 .and. 
     *   janei(k,1).eq.1 .and. jnnei(k,1).eq.3 .and. 
     *   janei(k,2).eq.16 .and. jnnei(k,2).eq.1 ) goto 6016 
c 
      if( nnei(k).eq.4 .and. 
     *   janei(k,1).eq.1 .and. jnnei(k,1).eq.2 .and. 
     *   janei(k,2).eq.6 .and. jnnei(k,2).eq.1 .and.
     *   janei(k,3).eq.16 .and. jnnei(k,2).eq.1 ) goto 6016 
c
c 
      if( nnei(k).eq.3 .and. 
     *   janei(k,1).eq.6 .and. jnnei(k,1).eq.1 .and. 
     *   janei(k,2).eq.7 .and. jnnei(k,2).eq.1 .and.
     *   janei(k,3).eq.16 .and. jnnei(k,2).eq.1 ) goto 6018 
cKNJ
      return
cKNJstart
6095  continue
      iexp(k) = 1
      write(nout,97) 'c2Cb(?)'
      ij=0
      do 6096 i=1,nnei(k)
        iat = jnei(k,i)
        ij = ij + 1
        if(ij.eq.1 .and. itype(iat).eq.5) ia1 = iat
        if(ij.eq.2 .and. itype(iat).eq.5) ia1 = iat
        if(ij.eq.3 .and. itype(iat).eq.5) ia1 = iat
6096   continue 
      if(itype(ia1).eq.5) then
        do i=1,nnei(ia1)
          jat = jnei(ia1,i)
          if( jat.ne.k ) then
            if( itype(jat).eq.9) then
              group(k) = 'C-b-F3'
              return
            else
              write(nout,999)  '(?)'
              call PrintGroup(k)
              goto 111
            endif
          endif
        enddo          
       endif
       goto 111
c
6070  continue 
      iexp(k) = 1
      write(nout,97) 'Ccno(?)'
      ij = 0
      do 6071 i=1,nnei(k)
        iat = jnei(k,i)
        ij = ij + 1
        if( ij.eq.1 ) ia1 = iat
        if( ij.eq.2 ) ia2 = iat
        if( ij.eq.3 ) ia3 = iat
c
c     checking which neighbour is oxygen and which is nitrogen
c
6071  continue 
c     we want to know if C and N are in the same planar ring
c
      if( itype(ia1).eq.7 .and. itype(ia2).eq.8) then
      iN = ia1
      iO = ia2
      call IsItCyclic(iN,lcyc,lPG,lPG2,iNr,mem)
      call IsItCyclic(iO,lcyc,lPG,lPG2,iOr,mem)       
      if( lcyc .and. lPG .and. (iCr.eq.iNr)) goto 6072
      else if(itype(ia1).eq.8 .and. itype(ia2).eq.7) then
      iO = ia1
      iN = ia2
      call IsItCyclic(iN,lcyc,lPG,lPG2,iNr,mem)
      call IsItCyclic(iO,lcyc,lPG,lPG2,iOr,mem)       
      if( lcyc .and. lPG .and. (iCr.eq.iNr)) goto 6072
      else
      write(nout,*) ':('
      goto 6039      
      endif  
c
6072  continue
c
c     C and N in the same ring, O out of the ring
c
      if( iCr.ne.iOr ) then
        if( nnei(iO).eq.1 ) then
        group(k) = 'amido'
        return
        endif 
        if( nnei(iO).eq.2 ) then
        group(k) = 'amidoH'
        return
        endif
      endif
c
c     C, N and O in the same planar ring 
c
      if ( iCr.eq.iOr ) then
      group(k) = 'o-Cc-n'
      return
      endif
c 
cKNJend
c   
5999  continue 
c      write(nout,*) nnei(k), lcyc, lPG
c
      return
c
c      o2               o2 -  
c     //               //    
c  c4-C1--o3-x   vs c4-C1--o3  
c     
6000  continue     
c     we have cCoo(?), now figure out if one of the oxygens 
c     has H or C,N.. atoms
      iexp(k) = 1
      write(nout,97) 'Coo(?)'
c     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 6001 i=1,nnei(k)
        iat = jnei(k,i)
        if( itype(iat).ne.8 ) goto 6001
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
6001  continue        
c      write(nout,*) ' i1O, i2O = ',i1O,i2O
c     maximum number of nei. around attached O-atoms
      maxn = max( nnei(i1O), nnei(i2O) )
c     carboxylate - each O-atom is bonded to only 1 atom, i.e. k-th Carbon
      if( maxn.eq.1 ) then
        group(k) = 'carboxylate'
        return
      endif
c     identify O-atom which has more neighbors  
      iat = 0      
      if( nnei(i1O).eq.2 ) iat = i1O
      if( nnei(i2O).eq.2 ) iat = i2O
c      write(nout,*) ' nnei(i1O) = ',nnei(i1O)
c      write(nout,*) ' nnei(i2O) = ',nnei(i2O)
c      write(nout,*) ' iat = ',iat
      if( iat.eq.0 ) then
        write(nout,999)  'cOO(?)'
        call PrintGroup(k)
        goto 111
      endif
c.....working with O-atom connected to 2 neighbors
c     1st neighbor - H-atom
      if( janei(iat,1).eq.1 ) then
        group(k) = 'carboxylic'
        return
      endif
c     two carbons around O-atom       
      if( jtnei(iat).eq.1 .and. janei(iat,1).eq.6 ) then
c       identify 'other' carbon
        do i=1,nnei(iat)
          jat = jnei(iat,i)
          if( jat.ne.k ) then
            if( nnei(jat).eq.3) then
              group(k) = 'anhydride'
              return
            else if( nnei(jat).eq.4) then
              group(k) = 'ester'
              return
            else
              write(nout,999)  'cOO(?)'
              call PrintGroup(k)
              goto 111
            endif
          endif
        enddo          
      else
        write(nout,999)  'cOO(?)'
        call PrintGroup(k)
      endif      
      goto 111
cKNJ
6016  continue
      iexp(k) = 1
      write(nout,97) 'C-s(?)'
      do i=1,nnei(k)
        iat = jnei(k,i)
       if( itype(iat).eq.16 ) then
        if( janei(iat,1).eq.30 .or.janei(iat,2).eq.30 ) then
!           do 6017 j=1, nnei(iat)
!           if( j.eq.1 ) ia1 = janei(iat,1)
!           if( j.eq.2 ) ia2 = janei(iat,2)
!        endif
!6017  continue 
!           if( ia1.eq.30 .or. ia2.eq.30) then
            group(k) = 'C-s-Zn'
            return
            else
              write(nout,999)  'c-s-(?)'
              call PrintGroup(k)
             goto 111
           endif
        endif
        enddo
      goto 111
c 
6018  continue
      iexp(k) = 1
      write(nout,97) 'c-C(=s)-n(?)'
      do i=1,nnei(k)
        iat = jnei(k,i)
       if( itype(iat).eq.16 ) then
        if( nnei(iat).eq.1 ) then
            group(k) = 'THIO'
            return
cPK 13_02_2018 
        else if( nnei(iat).eq.4 ) then
            group(k) = 'SX4'
            return
            else if( nnei(iat).eq.2) then
            group(k) = 'THIO-X'
            return
            else
              write(nout,999)  'c-C(=s)-n(?)'
              call PrintGroup(k)
             goto 111
           endif
        endif
       enddo
      goto 111
cKNJ
cPMD
cPMD      o2               o2 -  
cPMD     //               //    
cPMD  n4-C1--o3-x   vs n4-C1--o3  
cPMD     
6005  continue     
cPMD     we have nCoo(?), now figure out if one of the oxygens 
cPMD     has H or C,N.. atoms
      iexp(k) = 1
      write(nout,97) 'Coo(?)'
cPMD     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 6006 i=1,nnei(k)
        iat = jnei(k,i)
        if( itype(iat).ne.8 ) goto 6006
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
6006  continue        
cPMD      write(nout,*) ' i1O, i2O = ',i1O,i2O
cPMD     maximum number of nei. around attached O-atoms
      maxn = max( nnei(i1O), nnei(i2O) )
cPMD     carboxylate - each O-atom is bonded to only 1 atom, i.e. k-th Carbon
      if( maxn.eq.1 ) then
        group(k) = 'carbamate'
        return
      endif
cPMD     identify O-atom which has more neighbors  
      iat = 0      
      if( nnei(i1O).eq.2 ) iat = i1O
      if( nnei(i2O).eq.2 ) iat = i2O
cPMD      write(nout,*) ' nnei(i1O) = ',nnei(i1O)
cPMD      write(nout,*) ' nnei(i2O) = ',nnei(i2O)
cPMD      write(nout,*) ' iat = ',iat
      if( iat.eq.0 ) then
        write(nout,999)  'cOO(?)'
        call PrintGroup(k)
        goto 111
      endif
cPMD.....working with O-atom connected to 2 neighbors
cPMD     1st neighbor - H-atom
      if( janei(iat,1).eq.1 ) then
        group(k) = 'carbamic acid'
        return
      endif
cPMD     two carbons around O-atom       
      if( jtnei(iat).eq.1 .and. janei(iat,1).eq.6 ) then
cPMD       identify 'other' carbon
        do i=1,nnei(iat)
          jat = jnei(iat,i)
          if( jat.ne.k ) then
            if( nnei(jat).eq.4) then
              group(k) = 'carbamate ester'
              return
            else
              write(nout,999)  'cOO(?)'
              call PrintGroup(k)
              goto 111
            endif
          endif
        enddo          
      else
        write(nout,999)  'cOO(?)'
        call PrintGroup(k)
      endif      
      goto 111
cPMD end
c
c      c/h  
c       | (cyc)
c    n -C- c/h
c       |       
c      c/h
c
6010  continue
      iexp(k) = 1
      write(nout,97) 'n-C-(c/h)3'
c     check for "other" (i.e. NOT C and NOT H) atom types around carbon
      do i=1,jtnei(k)-1
        if( janei(k,i).ne.1 .and. janei(k,i).ne.6 ) then
          write(nout,999) 'n-C-(x/y/z)'
          goto 111
        endif
      enddo
c       
      iN = jsnei(k,jtnei(k),1)
      write(nout,96) ' Checking connected N-atom = ',atom(iN)
      if( nnei(iN).eq.4 ) then
        group(k) = 'ammonium'
        return
      endif
      if( nnei(iN).eq.3 ) then
c       check planarity of N-atom-group
        call CollectNei(iN,id,imem)
        call Planarity(0,lPG,imem,id,3)
        if( lPG ) then
          group(k) = 'planar amine'
        else
          group(k) = 'pyramidal amine'
        endif
        return        
      endif
      write(nout,999) 'n-C-(c/h)3'
      goto 111
c   ___
c  // \\
c       C-o-x
c  \___/
c   ---
c
6020  continue
      iexp(k) = 1
      write(nout,97) 'C(cyc)-o-x'
      iO =  jsnei(k,jtnei(k),1) 
      call IsItCyclic(iO,lcyc,lPG, lPG2,iOr,mem)
c     check if oxygen is in the same ring
      if( lcyc .and. lPG .and. (iCr.eq.iOr) ) goto 6029
c     check environment of oxygen
      if( nnei(iO).eq.1 ) then
        group(k) = 'phenolate'
        return
      endif
      if( nnei(iO).eq.2 ) then
        if( janei(iO,1).eq.1 .or. jtnei(iO).eq.1 ) then   
          group(k) = 'phenolic'
          return
        endif  
      endif 
      goto 111
cKNJstart
c     oxygen in the same ring
6029  continue
c     check environment of oxygen
      if( nnei(iO).eq.2 .and. janei(iO,1).eq.6 .and. 
     *  jtnei(iO).eq.1) then
        group(k) = 'furano'
        return
      endif
      if( nnei(iO).eq.2 .and. janei(iO,1).eq.6 .and. 
     *  janei(iO,2).eq.7) then
        group(k) = 'furanoN'
        return
      endif
      goto 111
cKNJ
6080  continue
      iexp(k) = 1
      write(nout,97) 'Cn2(cyc)=o/Cno(cyc)-n(?)'
      iO =  jsnei(k,jtnei(k),1) 
      call IsItCyclic(iO,lcyc,lPG,lPG2,iOr,mem)
c     check if oxygen is in the same ring
      if( lcyc .and. lPG .and. (iCr.eq.iOr) ) goto 6089
c     check environment of oxygen
      if( nnei(iO).eq.1 ) then
        group(k) = 'aromaturea'
        return
      endif
      goto 111
cKNJ
c     oxygen in the same ring
6089  continue
      if( nnei(iO).eq.2 ) then
        group(k) = 'pirofurano'
        return
      endif
      goto 111
cKNJend
c
c   ___              ___    
c  // \\   /        // \\   
c       C-n     or      C-h/c       ring is planar
c  \___/   \        \__n/  
c   ---
c
6030  continue
      iexp(k) = 1
      write(nout,97) 'C(cyc)-n='
      iN =  jsnei(k,jtnei(k),1) 
      call IsItCyclic(iN,lcyc,lPG,lPG2,iNr,mem)
cPMD  check if nitrogen is in the same ring and has 3 c/h neighbors
cPMD  nsp2(3)-=C(c/h)=-c
      if( lcyc .and. lPG .and. (iCr.eq.iNr) .and. nnei(iN).eq.3 ) then
        do j=1,nnei(iN)
          jat = jnei(iN,j)
          if( itype(jat).ne.1 .and. itype(jat).ne.6 ) goto 6039 
        enddo
        group(k) = 'substituted nitrogen'
        return
      endif      
cPMD  check if nitrogen is in the same ring and has 2 c as neighbors
cPMD  nsp2(2)-=C(c/h)=-c
      if( lcyc .and. lPG .and. (iCr.eq.iNr) .and. nnei(iN).eq.2 ) then
        group(k) = 'non-substituted nitrogen'
        return
      endif      
cPMD check if nitrogen is in the same ring (and has more than 3 neighbors)
      if( lcyc .and. lPG .and. (iCr.eq.iNr) ) goto 6039
cPMD end
c     check for =C-no2 group
      if( nnei(iN).eq.3 .and. janei(iN,jtnei(iN)).eq.8 .and.
     *                        jnnei(iN,jtnei(iN)).eq.2 ) then
        group(k) = 'nitrophenyl'
        return
      endif
c
cKNJ
c     check for =C-noh group
      if( nnei(iN).eq.2 .and. janei(iN,jtnei(iN)).eq.8 .and.
     *                        jnnei(iN,jtnei(iN)).eq.1 ) then
        group(k) = 'ph-noh'
        return
      endif
cKNJ
c     check for =C-n=(c/h) group
      if( nnei(iN).eq.3 ) then
        do j=1,nnei(iN)
          jat = jnei(iN,j)
          if( itype(jat).ne.1 .and. itype(jat).ne.6 ) goto 6039
        enddo
        group(k) = 'aminophenyl'
        return
      endif
cPMD  check for =C-n=(c/h)+ group
      if( nnei(iN).eq.4 ) then
        do j=1,nnei(iN)
          jat = jnei(iN,j)
          if( itype(jat).ne.1 .and. itype(jat).ne.6 ) goto 6039
        enddo
        group(k) = 'protonated aminophenyl'
        return
      endif
c
cPMD
c   ___n             
c  // \\           
c       C-x            ring is planar
c  \__n/           
c   ---
cKNJ  C-n= (C in ring, 3 n as a neighbor)
      if( nnei(k).eq.3 .and. lcyc .and. lPG  .and.       
     *   janei(k,jtnei(k)).eq.7 .and. jnnei(k,jtnei(k)).eq.3 ) goto 6060 
cKNJ
      goto 111
c
6031  continue 
      iexp(k) = 1
      write(nout,97) 'C(2cyc)-n='
      iN =  jsnei(k,jtnei(k),1) 
      call IsItCyclic(iN,lcyc,lPG,lPG2,iNr,mem)
      if( lcyc .and. lPG .and. mem.eq.5 ) then
        group(k) = 'C(2cyc)-n5'
        return
      endif      
      goto 111 
c
6035  continue
      iexp(k) = 1
      write(nout,97) 'n-C(cyc)-n'
      ij = 0
      do 6036 i=1,nnei(k)
        iat = jnei(k,i)
cPMD/cKNJstart  check if both nitrogen are in the same ring 
        call IsItCyclic(iat,lcyc,lPG,lPG2,iatr,mem)
      write(nout,*) ' iatr, iCr = ',iatr, iCr
        if( itype(iat).ne.7 ) goto 6036
        if( itype(iat).eq.7 .and. .not. ( lcyc .and. lPG )) goto 6090
        ij = ij + 1
        if( ij.eq.1 .and. lcyc .and. lPG .and. (iCr.eq.iatr) .and.
     *  itype(iat).eq.7 ) i1N = iat
        if( ij.eq.2 .and. lcyc .and. lPG .and. (iCr.eq.iatr) .and.
     *  itype(iat).eq.7 ) i2N = iat
        if( ij.eq.1 .and. lcyc .and. lPG .and. (iCr.ne.iatr) .and.
     *  itype(iat).eq.7 ) then
          in1 = jnei(k,1)
          in2 = jnei(k,2)
          in3 = jnei(k,3)
          nn = nnei(in1) + nnei(in2) + nnei(in3)
          if ( nn.eq.8 ) then
          group(k) = 'NNring1'
          else
          group(k) = 'NNring2'
          endif
        return
        goto 111
        endif
        if( ij.eq.2 .and. lcyc .and. lPG .and. (iCr.ne.iatr) .and.
     *  itype(iat).eq.7 ) then
          in1 = jnei(k,1)
          in2 = jnei(k,2)
          in3 = jnei(k,3)
          nn = nnei(in1) + nnei(in2) + nnei(in3)
          if ( nn.eq.8 ) then
          group(k) = 'NNring1'
          else
          group(k) = 'NNring2'
          endif
        return
        goto 111
        endif
c
cPMD  nsp2(3)-=C(c/h)=-nsp2(3)
6036  continue    
      write(nout,*) ' i1N, i2N = ',i1N,i2N
      write(nout,*) ' i1Nn, i2Nn = ',nnei(i1N), nnei(i2N)
cPMD     maximum number of nei. around attached N-atoms
         maxn = max( nnei(i1N), nnei(i2N) )
cKNJ
         if( maxn.eq.2 ) then
         group(k) = 'non-substituted nitrogens'
         return
         endif
c
         minn = min( nnei(i1N), nnei(i2N) )
c
         if( minn.eq.3 ) then
          if(janei(i1N,jtnei(i1N)).eq.30 ) then
              group(k) = 'ZNnCn'
              return
          endif
c
          if(janei(i2N,jtnei(i2N)).eq.30 ) then
              group(k) = 'ZNnCn'
              return
          endif
c
         group(k) = 'substituted nitrogens'
         return
         endif
cKNJ
         if( minn.eq.2 .and. maxn.eq.3 ) then
         group(k) = 'mixed nitrogens'         
         return
         endif
      goto 111
cKNJ
6090  continue
      write(nout,97) ' n-=Cn=-c (?) '
         ia1 = jnei(k,1)
         ia2 = jnei(k,2)
         ia3 = jnei(k,3)
         maxn = max( nnei(ia1), nnei(ia2), nnei(ia3) )
         minn = min( nnei(ia1), nnei(ia2), nnei(ia3) )
         sumn = 0
         sumn = sumn + nnei(ia1) + nnei(ia2) + nnei(ia3)         
cKNJ     check for the nitro group
         if( janei(ia1,jtnei(ia1)).eq.8 .or. janei(ia2,jtnei(ia2)).eq.8 
     *   .or. janei(ia3,jtnei(ia3)).eq.8) goto 6091
cKNJ     no oxygens attached
         if( janei(ia1,jtnei(ia1)).ne.8 .and. janei(ia2,jtnei(ia2)).ne.8 
     *   .and. janei(ia3,jtnei(ia3)).ne.8) goto 6092 
cKNJ
      goto 111
6091  continue
      write(nout,*) ' ia1, ia2, ia3 = ',ia1,ia2,ia3
cKNJ   
         if( maxn.eq.3 .and. minn.eq.2 .and. sumn.eq.8 ) then
           group(k) = 'nitro-n2n3c3'
         return
         endif
cKNJ
          if( maxn.eq.3 .and. minn.eq.3 .and. sumn.eq.9 ) then
          group(k) = 'nitro-n3n3c3'
         return
         endif
      goto 111
cKNJ
6092  continue  
      write(nout,*) ' ia4, ia5, ia6 = ',ia1,ia2,ia3
cKNJ   
         if( maxn.eq.3 .and. minn.eq.2 .and. sumn.eq.8 ) then
           group(k) = 'n2n3c3'
         return
         endif
cKNJ
         if( maxn.eq.3 .and. minn.eq.3 .and. sumn.eq.9 ) then
           group(k) = 'n3n3c3'
         return
         endif
      goto 111
cKNJ
cPMDcKNJend
6039  continue
      write(nout,999) 'C(cyc)-n='
      goto 111

c   ---Cany
c  /    \
c        C = x     ring is planar
c  \    /
c   ---Csp3
6040  continue
      iexp(k) = 1
      write(nout,97) 'c(sp3,cyc)-C(cyc)=x'
      do 6049 i=1,nnei(k)
        iat = jnei(k,i)
        if( itype(iat).ne.6 ) goto 6049
        call IsItCyclic(iat,lcyc,lPG,lPG2,iC1,mem)
        if( nnei(iat).eq.4 .and. iC1.eq.iCr ) then
          group(k) = 'cycloene'
          return
        endif
      goto 111
6049  continue
      write(nout,999) 'c(sp3,cyc)-C(cyc)=x'
c     this criterion is too loose, so we want to go back and check
c     another exceptions.... 
      iexp(k) = 0
      goto 5999
cPMD start
6050  continue
      iexp(k) = 1
      write(nout,97) 'c2Co'
      iO =  jsnei(k,jtnei(k),1) 
c     check environment of oxygen
      if( nnei(iO).eq.1 ) then
        group(k) = 'ketone'
        return
      else
        write(nout,999) 'c2--C--o-'
        call PrintGroup(k)
        goto 111
      endif
cKNJ
cKNJ  3 nitrogens as neighbours
6060  continue
      iexp(k) = 1
      write(nout,97) 'n-Cn-n'
      ij = 0
      do 6066 i=1,nnei(k)
        iat = jnei(k,i)
        call IsItCyclic(iat,lcyc,lPG, lPG2,iatr,mem)
        ij = ij + 1
        if( ij.eq.1 ) i1N = iat
        if( ij.eq.2 ) i2N = iat
        if( ij.eq.3 ) i3N = iat
6066  continue    
cKNJ     write(nout,*) ' i1N, i2N = ',i1N,i2N
cKNJ     maximum number of nei. around attached N-atoms
         maxn = max( nnei(i1N), nnei(i2N), nnei(i3N) )
         minn = min( nnei(i1N), nnei(i2N), nnei(i3N) )
         sumn = 0
         sumn = sumn + nnei(i1N) + nnei(i2N) + nnei(i3N)         
c
         if( maxn.eq.3 .and. sumn.eq.7 ) then
         group(k) = '3N-sp2(232)'
         return
         endif
c
         if( minn.eq.3 ) then
         group(k) = '3N-sp2(333)'
         return
         endif
c        
         if( maxn.eq.3 .and. sumn.eq.8 ) then 
         group(k) = '3N-sp2(233)'         
         return
         endif
cPMDcKNJend
cPMD start
c---- NITROGEN atom ------------------------------------------------------
c
7     continue
c     c-Nc1-c1 N-planar : c1=csp3 or csp2 ?
      call CollectNei(k,id,imem)
      lPG = .false.
      if(imem.ge.3) call Planarity(0,lPG,imem,id,3)
      if( lPG .and. nnei(k).eq.3 .and. jtnei(k).eq.1 
     *                           .and. janei(k,1).eq.6 ) goto 7010
cKNJ
c     N-O : c-N=o or c-N-o-h group ?
      if( nnei(k).eq.2 .and. janei(k,1).eq.6
     *                 .and. janei(k,2).eq.8 ) goto 7020
cKNJ
      return
c
7010  continue
c     we have c-Nc1-c1 N-planar, now figure out the valency of carbons
      iexp(k) = 1
      write(nout,97) 'c-Nc1-c1 N-planar c1=?'
c     figure out the seq. numbers of connected C-atoms
      ij = 0
      do 7011 i=1,nnei(k)
        iat = jnei(k,i)
        if( itype(iat).ne.6 ) goto 7011
        ij = ij + 1
        if( ij.eq.1 ) i1C = iat
        if( ij.eq.2 ) i2C = iat
	if( ij.eq.3 ) i3C = iat
7011  continue        
      sumn = 0
      sumn = sumn + nnei(i1C) + nnei(i2C) + nnei(i3C)
      minN = min( nnei(i1C), nnei(i2C), nnei(i3C) )
      maxN = max( nnei(i1C), nnei(i2C), nnei(i3C) )
c     c-Nc1-c1 N-planar : c=csp2 c1=csp3
      if( sumn.eq.11 .and. minN.eq.3 .and. maxN.eq.4 ) then
        group(k) = 'csp2sp3sp3'
        return
c     c-Nc1-c1 N-planar : c=csp3 c1=csp2
      else if( sumn.eq.10 .and. minN.eq.3 .and. maxN.eq.4 ) then
        group(k) = 'csp3sp2sp2'
        return
      else
        write(nout,999)  'c-Nc1-c1 planar'
        call PrintGroup(k)
        goto 111
      endif
cPMD end     
cKNJ
7020  continue
c     we have c--N--o
      iexp(k)=1
      write(nout,97) 'c--N--o(?)'
c     get environment of O-atom
      iO =  jsnei(k,jtnei(k),1)       
      if( itype(iO).eq.8 ) then 
        if( nnei(iO).eq.2 .and. janei(iO,1).eq.1
     *                    .and. janei(iO,2).eq.7 ) then         
        group(k) = 'N-OH'
        return
cKNJ        
        else if( nnei(iO).eq.1 .and. janei(iO,1).eq.7 ) then       
        group(k) = 'N=O'
        return
cKNJ        
        else
        write(nout,999)  'c--N--c?'
        call PrintGroup(k)
        goto 111
        endif
cKNJ
      endif
cKNJ
c
c---- OXYGEN atom ------------------------------------------------------
c
8     continue
c     c-O-h : alcohol/phenol or carboxyl group ?
      if( nnei(k).eq.2 .and. janei(k,1).eq.1 
     *                 .and. janei(k,2).eq.6 ) goto 8010
cKNJ     O=c 
      if( nnei(k).eq.1 .and. janei(k,1).eq.6 ) goto 8020
c     c-O-c : anhydride or ester or carbamate ester
      if( nnei(k).eq.2 .and. jtnei(k).eq.1 
     *                 .and. janei(k,1).eq.6 ) goto 8030
cKNJ  n-O : nitro1 or nitro2 group ?
      if( nnei(k).eq.1 .and. janei(k,1).eq.7 ) goto 8060
cKNJ
cPMD  O=p : in cpo3 or cpo2(o-x) or cpo(o-x)2 or 
cKNJ: c3pO or pO4conf
      if( nnei(k).eq.1 .and. janei(k,1).eq.15 ) goto 8040
cPMD
      return
c
c     ***  ?c-O-h  : ALCOHOL or CARBOXYLIC group ?  ***
c
8010  iexp(k)=1
      write(nout,97) 'c-O-h'
c     what is x-atom :    
      iC = jsnei(k,2,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
c     serial number of oxygen atom type among neighbors of carbon 
      jox = 0
      do i=1,jtnei(iC)      
        if( janei(iC,i).eq.8 ) jox = i
      enddo
c     jot = number of o-atoms around of c-atom
      jot = jnnei(iC,jox)
c     only 1 oxygen:
      if( jot.eq.1 ) then
c       c-atom is 3-coordinated and in a planar ring
        if( nnei(iC).eq.3 .and. lcyc .and. lPG ) then
          group(k) = 'phenol'     
          return
        endif  
c       c-atom is 4-coordinated or 3-coordinated but not in the ring
        if( nnei(iC).eq.3 .or. nnei(iC).eq.4 ) then
          group(k) = 'alcohol'     
          return
        endif
        write(nout,999)  'c-O-h'
        call PrintGroup(k)
        goto 111
      endif
c     2 oxygens around carbon
c     must check how many bonds the "other" oxygen has...
c     first identify the "other" oxygen     
      if( jot.eq.2 ) then       
        iO = 0
        do j=1,jnnei(iC,jox)
          ja = jsnei(iC,jox,j)
          if( ja.ne.k ) iO = ja 
        enddo
        if( iO.eq.0 ) call Err('DetExcept',
     *              'Could not identify "other" O-atom')
c       if "other" O-atom has only 1 neighbor, 
c       then THIS oxygen is CARBOXYLIC
        if( nnei(iO).eq.1 ) then
          group(k) = 'carboxylic'  
          return 
        endif
c       otherwise - ALCOHOL 
        group(k) = 'alcohol'  
        return   
      endif   
c    "crazy" number of o-atoms around carbon       
      write(nout,999) 'x-O-h'
      call PrintGroup(k)
      goto 111
c
c     *** O=c  : aldehyde, carboxyl, carboxylate, carbonyl, ketone ***
c     ***        ester, anhydride...                               ***
cKNJ
8020  iexp(k)=1
      write(nout,97) '(?)c=O'
c     get environment of c-atom 
      iC = jnei(k,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
      if( nnei(iC).ne.3 ) then
         write(nout,999) '(?)c=O [1]' 
         call PrintGroup(k)
         goto 111
      endif   
999   format(/' DetExcept :: Unknown ',a,' functional group',//,
     *        ' See file group.ang for more info...'/)   
c     quasiAROMATIC: c-atom is cyclic and planar 
c     serial number of oxygen atom type among neighbors of carbon 
c       c-atom is 3-coordinated and in a planar ring
      if( nnei(iC).eq.3 .and. lcyc .and. lPG ) then
        if( janei(iC,1).eq.7 .and. jnnei(iC,1).eq.2 ) then
        group(k)='aromaturea'
        return
        endif
c
        if( janei(iC,1).eq.6 .and. jnnei(iC,1).eq.2 ) then
        group(k)='aromatketone'
        return
        endif
c
        if( janei(iC,1).eq.6 .and. jnnei(iC,1).eq.1 .and.
     *    janei(iC,2).eq.8 .and. jnnei(iC,2).eq.2 ) then
          group(k) = 'aromatester'     
          return
        endif  
c
        if( janei(iC,1).eq.6 .and. janei(iC,2).eq.7 .and. 
     *    janei(iC,3).eq.8 ) then
          group(k)='aromatamide'
        return
        endif         
       endif
cKNJ
c     ALDEHYDE : c-atom has to be bonded to c- and h-atoms    
c
      if( janei(iC,1).eq.1 .and. janei(iC,2).eq.6 ) then
        group(k)='aldehyde'
        return
      endif     
c
c     KETONE : c-atom is bonded to 2 other c-atoms
c
      if( janei(iC,1).eq.6 .and. jnnei(iC,1).eq.2 ) then
        group(k)='ketone'
        return
      endif     
c
cKNJ  UREA : c-atom is bonded to 2 nitrogen atoms
c
      if( janei(iC,1).eq.7 .and. jnnei(iC,1).eq.2 ) then
        group(k)='urea'
        return
      endif
cKNJ      
c                O1    o5           O1                              O1    
c     ANHYDRIDE  ||    ||  or ESTER ||     |   or CARBOXYL(ATE/IC)  || 
c              c-c2-o3-c4-        c-c2-o3-c4-                     c-c2-o3-(h)
c                                         |
c
      if( janei(iC,1).eq.6 .and. jnnei(iC,1).eq.1 .and.
     *    janei(iC,2).eq.8 .and. jnnei(iC,2).eq.2 ) then
c       find O3
        iO3 = 0
        do i=1,jnnei(iC,2)
          if( jsnei(iC,2,i) .ne. k ) iO3 = jsnei(iC,2,i)
        enddo
c       if O3 has only 1 neighbor - carboxylate        
        if( nnei(iO3).eq.1) then
          group(k)='carboxylate'
          return
        endif
c       if O3 has 2 neigbors and 1st one is H-atom - carboxylic
        if( nnei(iO3).eq.2 .and. janei(iO3,1).eq.1 ) then
          group(k)='carboxylic'
          return
        endif
c        if( nnei(iO3).ne.2 .or. jtnei(iO3).ne.1 .or. 
c     *     janei(iO3,1).ne.6 ) goto 111
c        write(nout,*) ' iO3 = ',iO3
c       find C4
        iC4 = 0
        do i=1,nnei(iO3)
          if( jnei(iO3,i).ne.iC ) iC4 = jnei(iO3,i)
        enddo          
c        write(nout,*) ' iC4 = ',iC4
c      if C4 is 3-coordinated and has 2 oxygens - anhydride
       if( nnei(iC4).eq.3 ) then
         ino = 0
         do i=1,jtnei(iC4)
           if( janei(iC4,i).eq.8 ) ino = jnnei(iC4,i)
         enddo         
         if( ino.ne.2 ) goto 111
         group(k) = 'anhydride'
         return
       endif
c      if C4 is 4-coordinated - ester
       if( nnei(iC4).eq.4 ) then
         group(k) = 'ester'
         return         
       endif
c      Unknown  c=O group 
       write(nout,999) '(?)c=O [2]' 
       call PrintGroup(k)
       goto 111       
      endif
cPMD                 O1                                        O1    
cPMD CARBAMATE ESTER ||     |   or  CARBAMIC ACID/CARBAMATE    || 
cPMD               n-c2-o3-c4-                               n-c2-o3-(h)
cPMD                        |
cPMD
      if( janei(iC,1).eq.7 .and. jnnei(iC,1).eq.1 .and.
     *    janei(iC,2).eq.8 .and. jnnei(iC,2).eq.2 ) then
cPMD       find O3
        iO3 = 0
        do i=1,jnnei(iC,2)
          if( jsnei(iC,2,i) .ne. k ) iO3 = jsnei(iC,2,i)
        enddo
cPMD       if O3 has only 1 neighbor - carbamate        
        if( nnei(iO3).eq.1) then
          group(k)='carbamate'
          return
        endif
cPMD       if O3 has 2 neigbors and 1st one is H-atom - carbamic acid
        if( nnei(iO3).eq.2 .and. janei(iO3,1).eq.1 ) then
          group(k)='carbamic acid'
          return
        endif
cPMD        if( nnei(iO3).ne.2 .or. jtnei(iO3).ne.1 .or. 
cPMD     *     janei(iO3,1).ne.6 ) goto 111
cPMD        write(nout,*) ' iO3 = ',iO3
cPMD       find C4
        iC4 = 0
        do i=1,nnei(iO3)
          if( jnei(iO3,i).ne.iC ) iC4 = jnei(iO3,i)
        enddo          
cPMD        write(nout,*) ' iC4 = ',iC4
cPMD      if C4 is 4-coordinated - carbamate ester
       if( nnei(iC4).eq.4 ) then
         group(k) = 'carbamate ester'
         return         
       endif
cPMD      Unknown  c=O group 
       write(nout,999) '(?)c=O [2]' 
       call PrintGroup(k)
       goto 111       
      endif
cPMD end
c
c     THIOESTER : r-s-cO-r'
c
      jCn = jtnei(iC)
      if( janei(iC,jCn).eq.16 .and. jnnei(iC,jCn).eq.1 ) then
        group(k)='thioester'
        return
      endif
c
c     AMIDE
c
      if( janei(iC,1).eq.6 .and. janei(iC,2).eq.7 .and. 
     *    janei(iC,3).eq.8 ) then
        iN = jsnei(iC,2,1)
        if( nnei(iN).ge.3 ) then
          group(k)='amide'
          return
        endif  
      endif
c
c     Unknown  c=O group 
c
      write(nout,999) '(?)c=O [3]' 
      call PrintGroup(k)
      goto 111
c
c     *** c-O-c : ESTER or ANHYDRIDE or ... ***
c
8030  iexp(k)=1
      write(nout,97) 'c-O-c'
      iC1 = jnei(k,1)
      iC2 = jnei(k,2)
      iC1nei = nnei(iC1)
      iC2nei = nnei(iC2)
c     checking for ANHYDRIDE : each c-atom has to be connected to 2 oxygens and 1 carbon
      if( iC1nei.eq.3 
     *          .and. janei(iC1,1).eq.6 .and. jnnei(iC1,1).eq.1
     *          .and. janei(iC1,2).eq.8 .and. jnnei(iC1,2).eq.2 ) then
         if( iC2nei.eq.3 
     *         .and. janei(iC2,1).eq.6 .and. jnnei(iC2,1).eq.1
     *         .and. janei(iC2,2).eq.8 .and. jnnei(iC2,2).eq.2 ) then
           group(k)='anhydride'
         return
         endif
      endif 
c     checking for ESTER and CARBAMATE ESTER
c     checking for ESTER
      maxN = max( iC1nei, iC2nei )
      minN = min( iC1nei, iC2nei )
c      write(nout,*) ' minN, maxN = ',minN, maxN
      if( minN.eq.3 .and. (maxN.eq.3.or.maxN.eq.4) ) then
        iflag = 0
c       loop over atom types around carbons
        do i=1,2
          iC = jnei(k,i)
c         in ESTER one of the c-atoms should be 3-coordinated and bound
c         to 2 oxygens and 1 carbon  
c          write(nout,*) i, iC, nnei(iC), jtnei(iC)           
          if( nnei(iC).eq.3 .and. jtnei(iC).eq.2   
     *       .and. janei(iC,1).eq.6 .and. jnnei(iC,1).eq.1 
     *       .and. janei(iC,2).eq.8 .and. jnnei(iC,2).eq.2 )  iflag = 1
        enddo  
        if( iflag.eq.1 ) then
          group(k)='ester'
          return
        endif  
      endif  
cPMD     checking for CARBAMIDE ESTER
      maxN = max( iC1nei, iC2nei )
      minN = min( iC1nei, iC2nei )
cPMD      write(nout,*) ' minN, maxN = ',minN, maxN
      if( minN.eq.3 .and. (maxN.eq.3.or.maxN.eq.4) ) then
        iflag = 0
cPMD       loop over atom types around carbons
        do i=1,2
          iC = jnei(k,i)
cPMD         in CARBAMATE ESTER one of the c-atoms should be 3-coordinated and bound
cPMD         to 2 oxygens and 1 nitrogen  
cPMD          write(nout,*) i, iC, nnei(iC), jtnei(iC)           
          if( nnei(iC).eq.3 .and. jtnei(iC).eq.2   
     *       .and. janei(iC,1).eq.7 .and. jnnei(iC,1).eq.1 
     *       .and. janei(iC,2).eq.8 .and. jnnei(iC,2).eq.2 )  iflag = 1
        enddo  
        if( iflag.eq.1 ) then
          group(k)='carbamate ester'
          return
        endif  
      endif  
cPMD end      
c     who knows what it is... 
      write(nout,999) 'c-O-c' 
      call PrintGroup(k)
      goto 111
c
cPMD start
8040  iexp(k)=1
      write(nout,97) '(?)p=O'
cPMD    get environment of p-atom 
      iP = jnei(k,1)
      if( nnei(iP).ne.4 ) then
         write(nout,998) '(?)p=O [1]' 
         call PrintGroup(k)
         goto 111
      endif   
998   format(/' DetExcept :: Unknown ',a,' functional group',//,
     *        ' See file group.ang for more info...'/)   
c
cKNJ   get environment of p-atom
      if( nnei(iP).eq.4 .and.
     *   janei(iP,1).eq.6 .and. jnnei(iP,1).eq.3 .and.
     *   janei(iP,2).eq.8 .and. jnnei(iP,2).eq.1 ) then
        group(k)='C3PO'
        return
      endif
cKNJ
      if( nnei(iP).eq.4 .and.
     * janei(iP,1).eq.8 .and. jnnei(iP,1).eq.4 ) goto 8050
cKNJ
      if( nnei(iP).eq.4 .and.
     * janei(iP,1).eq.6 .and. jnnei(iP,1).eq.1 .and.
     * janei(iP,2).eq.8 .and. jnnei(iP,2).eq.3 ) goto 8041
cKNJ
8050  continue                  
c     we have pOooo(?), now figure out the valency of oxygens
      iexp(iP) = 1
      write(nout,97) 'pOooo(?)'
c     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 8051 i=1,nnei(iP)
        iat = jnei(iP,i)
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
        if( ij.eq.3 ) i3O = iat
        if( ij.eq.4 ) i4O = iat
8051  continue        
      sumn = 0
      sumn = sumn + nnei(i1O) + nnei(i2O) + nnei(i3O) + nnei(i4O)
c     O-Po3
      if( sumn.eq.4 ) then
        group(k) = 'pO4'
        return
c     O-P(OX)3
      else if( sumn.eq.7 ) then
        group(k) = 'O-p(ox)3'
        return
c     O3-P(OX)
      else if( sumn.eq.5 ) then
        group(k) = 'O3-p(ox)'
        return
      else
        write(nout,999)  'pOooo(?)'
        call PrintGroup(k)
        goto 111
      endif      
cKNJ
c
c
c       O2                 O2                O2
c      //                 //                //
c  c5-p1=O4-    vs    c5-p1-o4-x   vs   c5-p1-o4-x
c     \\                 \\                \ 
c      O3-                O3-               o3-x
c                  
8041  continue                  
c     we have cpOoo(?), now figure out the valency of oxygens
      iexp(iP) = 1
      write(nout,97) 'cpOoo(?)'
c     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 8042 i=1,nnei(iP)
        iat = jnei(iP,i)
        if( itype(iat).ne.8 ) goto 8042
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
        if( ij.eq.3 ) i3O = iat
8042  continue        
      sumn = 0
      sumn = sumn + nnei(i1O) + nnei(i2O) + nnei(i3O)
c     cPo3h
      if( sumn.eq.4 ) then
        group(k) = 'CPO2(O-X)'
        return
c     cPo3h2
      else if( sumn.eq.5 ) then
        group(k) = 'CPO(O-X)2'
        return
c     cPo3 
      else if( sumn.eq.3 ) then
        group(k) = 'CPO3'
        return
      else
        write(nout,999)  'cpOoo(?)'
        call PrintGroup(k)
        goto 111
      endif
cKNJ
8060  continue
      iexp(k)=1
      write(nout,97) '(?)N-O'
c     get environment of N-atom 
      iN = jnei(k,1)
      if( nnei(iN).ne.2 .and. nnei(iN).ne.3 ) then
         write(nout,999) '(?)N-O [1]' 
         call PrintGroup(k)
         goto 111
      endif   
cKNJ
      if( nnei(iN).eq.2 .and. 
     *    janei(iN,1).eq.6 .and. jnnei(iN,1).eq.1 .and.
     *    janei(iN,2).eq.8 .and. jnnei(iN,2).eq.1 ) then
        group(k) = 'NO'
        return
      endif  
cKNJ     
      if( nnei(iN).eq.3 .and. 
     *    janei(iN,1).eq.6 .and. jnnei(iN,1).eq.1 .and.
     *    janei(iN,2).eq.8 .and. jnnei(iN,2).eq.2 ) then
        group(k) = 'NO2'
        return
      endif
cPK NO2N Exception
      if( nnei(iN).eq.3 .and. 
     *    janei(iN,1).eq.7 .and. jnnei(iN,1).eq.1 .and.
     *    janei(iN,2).eq.8 .and. jnnei(iN,2).eq.2 ) then
        group(k) = 'NO2N'
        return
      endif
cKNJ
c
cPMD end
c
cPMD start
cPMD---- PHOSPHORUS atom ----------------------------------------------------
c
15    continue
      if( nnei(k).eq.4 .and. 
     * janei(k,1).eq.6 .and. jnnei(k,1).eq.1 .and.
     * janei(k,2).eq.8 .and. jnnei(k,3).eq.3 ) goto 1500
cKNJ
      if( nnei(k).eq.4 .and. 
     * janei(k,1).eq.8 .and. jnnei(k,1).eq.4 ) goto 1510
cKNJ
      return
c       o2                 o2                o2
c      //                 //                //
c  c5-P1=o4-    vs    c5-P1-o4-x   vs   c5-P1-o4-x
c     \\                 \\                \ 
c      o3-                o3-               o3-x
c                  
1500  continue                  
c     we have cPooo(?), now figure out the valency of oxygens
      iexp(k) = 1
      write(nout,97) 'cPooo(?)'
c     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 1501 i=1,nnei(k)
        iat = jnei(k,i)
        if( itype(iat).ne.8 ) goto 1501
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
	if( ij.eq.3 ) i3O = iat
1501  continue        
      sumn = 0
      sumn = sumn + nnei(i1O) + nnei(i2O) + nnei(i3O)
c     cPo3h
      if( sumn.eq.4 ) then
        group(k) = 'CPO2(O-X)'
        return
c     cPo3h2
      else if( sumn.eq.5 ) then
        group(k) = 'CPO(O-X)2'
        return
c     cPo3 
      else if( sumn.eq.3 ) then
        group(k) = 'CPO3'
        return
      else
        write(nout,999)  'cPooo(?)'
        call PrintGroup(k)
        goto 111
      endif
c
cPMD end
c
cKNJ
1510  continue                  
c     we have Poooo(?), now figure out the valency of oxygens
      iexp(k) = 1
      write(nout,97) 'Poooo(?)'
c     figure out the seq. numbers of connected O-atoms
      ij = 0
      do 1511 i=1,nnei(k)
        iat = jnei(k,i)
        ij = ij + 1
        if( ij.eq.1 ) i1O = iat
        if( ij.eq.2 ) i2O = iat
        if( ij.eq.3 ) i3O = iat
        if( ij.eq.4 ) i4O = iat
1511  continue        
      sumn = 0
      sumn = sumn + nnei(i1O) + nnei(i2O) + nnei(i3O) + nnei(i4O)
c     o3-P(ox)
      if( sumn.eq.5 ) then
        group(k) = 'o3-P(ox)'
        return 
c     o3-P(ox)  
      else if( sumn.eq.6 ) then
        group(k) = 'o2-P(ox)2'
        return
c     o3-P(ox)  
      else if( sumn.eq.7 ) then
        group(k) = 'o-P(ox)3'      
        return
      else
        write(nout,999)  'Poooo(?)'
        call PrintGroup(k)
        goto 111
      endif
c
cKNJ---- SULPHUR atom -----------------------------------------------------
c
16    continue 
      if( nnei(k).eq.1 .and. janei(k,1).eq.6 ) goto 1600
c
1600  continue
      iexp(k)=1
      write(nout,97) '(?)c=S'
c     get environment of c-atom 
      iC = jnei(k,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
      if( nnei(iC).ne.3 ) then
         write(nout,999) '(?)c=S [1]' 
         call PrintGroup(k)
         goto 111
      endif   
      if( nnei(iC).eq.3 .and. lcyc .and. lPG ) then
        if( janei(iC,1).eq.7 .and. jnnei(iC,1).eq.2 ) then
        group(k)='PHthiourea'
        return
        endif
c
        if( janei(iC,1).eq.6 .and. jnnei(iC,1).eq.2 ) then
        group(k)='PH=S'
        return
        endif
c
        if( janei(iC,1).eq.6 .and. janei(iC,2).eq.7 .and. 
     *    janei(iC,3).eq.16 ) then
          group(k)='aromathioamide'
        return
        endif         
      endif
c
       goto 111
c
c
cKNJ---- FLUORINE atom ----------------------------------------------------
c
9     continue
      if( nnei(k).eq.1 .and. janei(k,1).eq.6 ) then
      iexp(k)=1
      write(nout,97) '(?)c-F'
c     get environment of c-atom 
      iC = jnei(k,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
       if( lcyc .and. lPG ) then
        group(k) = 'PH-F'
       return
       goto 111
       endif
      endif
c
      goto 111      
c
cKNJ---- CHLORINE atom ----------------------------------------------------
c
17    continue
      if( nnei(k).eq.1 .and. janei(k,1).eq.6 ) then
      iexp(k)=1
      write(nout,97) '(?)c-Cl'
c     get environment of c-atom 
      iC = jnei(k,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
       if( lcyc .and. lPG  ) then
        group(k) = 'PH-Cl'
       return
       goto 111
       endif
      endif
c
c     if group is not identified - remove exception
c
cKNJ---- BORINE atom ----------------------------------------------------
c
35    continue
      if( nnei(k).eq.1 .and. janei(k,1).eq.6 ) then
      iexp(k)=1
      write(nout,97) '(?)c-Br'
c     get environment of c-atom 
      iC = jnei(k,1)
      call IsItCyclic(iC,lcyc,lPG,lPG2,iCr,mem)
       if( lcyc .and. lPG  ) then
        group(k) = 'PH-Br'
       return
       goto 111
       endif
      endif
c
c     if group is not identified - remove exception
c
111   continue
      iexp(k) = 0
      group(k) = ' '
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prints out to file 'group.ang' neigbors of k-th atom
      subroutine PrintGroup(k)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      parameter( nda = 100) 
      dimension idum(nda)
      call izero(idum,nda)
      ip = 89
      open(ip,file='group.ang',form='formatted',status='unknown')
      write(ip,'(a)') 'ANGSTROM'
      nd = 1
      idum(1) = k
c     loop over neighbors of k-th atom
      do 1 i=1,nnei(k)
        iat = jnei(k,i)
c       loop over neigbors of i-th neigbor of k-th atom
        do 2 j=1,nnei(iat)
          jat = jnei(iat,j)
c         loop over neigbors of j-th neigbor of i-th atom           
          do 3 l=1,nnei(jat)
            lat = jnei(jat,l)
            call CheckAt(nda,nd,idum,lat,*3)
3         continue
          call CheckAt(nda,nd,idum,jat,*2)
2       continue
        call CheckAt(nda,nd,idum,iat,*1)
1     continue
c     write atoms to platon file
      do i=1,nd
        iat = idum(i)
        write(ip,45) atom(iat), (xc(iat,l),l=1,3)
      enddo  
45    format('ATOM ',a,2x,3f12.4)      
      close(ip)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CheckAt(nda,nd,idum,l,*)
      dimension idum(nda)
      do i=1,nd
        if( l.eq.idum(i) ) return 1
      enddo
      nd = nd + 1
      idum(nd) = l
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     determine cartesian coordinates of aux. atoms which define
c     the local coordinate system for k-th atom and write them to 
c     PLATON file
c     cartesian coordinates of origin  -> xcao
c                             1st atom -> xca1
c                             2nd atom -> xca2                          
      subroutine Loc2Plat(k,xca1,xca2,ax1,ax2)
      include 'lsdb.inc'
      dimension xcao(3), xca1(3), xca2(3)
      character axj*1, ax1, ax2
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      character ax(3), adum*10
      dimension rx1(3), rx2(3), rx3(3), iax(3)
      dimension xa(3,3)
      dimension xd(3), yd(3), zd(3)
      common /iplat/ iplat
      common /files/ imas, ilog, ixp, idb, icon
      common /ifra/ ifra
c      character dlab*6
c      common /dummy/ ndum, dlab(ndumax), xdum(ndumax,3)
      data ax /'X','Y','Z'/
c
      do i=1,3
       xcao(i) = xc(k,i)
       iax(i)=0
      enddo
c      ilog=47
c      write(ilog,348) k, atom(k)(:ils(atom(k)))
c348   format(//' Setting up local axes in PLATON for atom ',i6,
c     *         2x,a,' :')
      write(iplat,503) atom(k), xcao
c
      do j=1,3
        rx1(j) = xca1(j) - xcao(j)
        rx2(j) = xca2(j) - xcao(j)
      enddo
c.....third axis axis: Lets find new 3rd - it is [ 1st X 2nd ] :
      call VMult(rx1,rx2,rx3)
c      write(ilog,*) ' '      
c      write(ilog,328) ' rx1 = ',rx1,rlength(rx1)      
c      write(ilog,328) ' rx2 = ',rx2,rlength(rx2)      
c      write(ilog,328) ' rx3 = ',rx3,rlength(rx3)     
c.....second axis again
      call VMult(rx3,rx1,rx2)
c      write(ilog,*) ' '      
c      write(ilog,328) ' rx1 = ',rx1,rlength(rx1)      
c      write(ilog,328) ' rx2 = ',rx2,rlength(rx2)     
c      write(ilog,328) ' rx3 = ',rx3,rlength(rx3)      
328   format(a,3f9.4,3x,f9.4)
c.....now assign axes to the correct arrays
      do i=1,2
        if(i.eq.1) axj=ax1
        if(i.eq.2) axj=ax2 
        do j=1,3
          if(axj.eq.ax(j)) then
            iax(j)=1
            do l=1,3
              if(i.eq.1) xa(j,l)=rx1(l)
              if(i.eq.2) xa(j,l)=rx2(l)
            enddo 
          endif
        enddo
      enddo
      do i=1,3
        if(iax(i).eq.0) then
          do j=1,3
            xa(i,j)=rx3(j)
          enddo
        endif
      enddo       
      do j=1,3
        rx1(j)=xa(1,j)
        rx2(j)=xa(2,j)
        rx3(j)=xa(3,j)
      enddo
c---  PLATON check
      call GetDummy(xcao,rx1,rx2,rx3,xd,yd,zd)
c      call get_dummy_label(k,adum)
      call GetDummyLabel(k,atom(k),adum)
      is=ils(adum)
c      write(ilog,'(/a)') 
c     * ' Check for PLATON: cartesian coordinates of axes:'      
c      write( ilog,503) 'X'//adum(1:is), xd
c      write( ilog,503) 'Y'//adum(1:is), yd
c      write( ilog,503) 'Z'//adum(1:is), zd
      write(iplat,503) 'X'//adum(1:is), xd
      write(iplat,503) 'Y'//adum(1:is), yd
      write(iplat,503) 'Z'//adum(1:is), zd
503   format('ATOM ',a8,3x,3f12.4)      
98    format(a,4f10.4)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GetDummyLabel(k,atom,adum)
      character adum*(*), atom*(*), line*255
      line=' '
      adum=' '
c      ib=1
c      if(atom(2:2).ne.'(') ib=2
c      write(line,'(i6,a)') k, atom(1:ib)
c 010603: only number
      write(line,'(i6)') k
      ik=0
      do i=1,100
        if(line(i:i).ne.' ') then
          ik=ik+1
          adum(ik:ik)=line(i:i)
        endif
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get dummy atoms for PLATON to draw local coordinate systems
c     coordinates of dummy atoms rescaled to 'dn' ANG from parent atom
      subroutine GetDummy(xo,xn,yn,zn,xd,yd,zd)
      dimension xo(3),xn(3),yn(3),zn(3),xd(3),yd(3),zd(3)
      dimension xv(3)
      common /files/ imas, ilog, ixp, idb, icon 
c.....
      call rzero(xv,3)
c.....rescale to 0.2 Angstroms
      dn=0.2D0
c.....x
      d1=dn/rlength(xn)
      d2=dn/rlength(yn)
      d3=dn/rlength(zn)
      do j=1,3
        xd(j)=d1*xn(j)+xo(j)
        yd(j)=d2*yn(j)+xo(j)
        zd(j)=d3*zn(j)+xo(j)
        xv(j) = xv(j) + d1*xn(j) + d2*yn(j) + d3*zn(j)
      enddo       
c      write(ilog,'(/a,3f11.5)')  ' Check vector: ',xv    
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     vector multiplication x*y=z
      subroutine VMult(x,y,z)
      dimension x(3),y(3),z(3)
      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=x(3)*y(1)-x(1)*y(3)
      z(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
c
c***********************************************************************
c
c                   D  A  T  A  B  A  N  K
c
c***********************************************************************
c
c read COPPENS LAB DATABANK and find pseudoatom parameters
c if there is no match, then assign default parameters :
c    Pval = free atom, kappas = 1, Plms = 0.0
c
      subroutine ReadDB(ik)
      include 'lsdb.inc'
      character atom(nat)*8, asym(nat)*5, asymr*5, asymi*5, groupj*40
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /files/ imas, io, ixp, idb, icon
      common /verbosity/ iverbose
c.....neighbor info from FindNei
      dimension iatn(maxnei), rnd(maxnei) 
      dimension ndt(maxnei), natnum(maxnei), ndtnum(maxnei,maxnei)
      common /neig/ inei, iatn, rnd, ity, ndt, natnum, ndtnum
c
      common /forH/ isp(nat), ihp(nat), ihnei(nat)
      common /forRing/ iring, imemb
c
c these parameters are necessary to make search - they can be input
c in this subroutine by other means
c  inei  - number of neighbours
c  iatn  - sequence numbers of neighbors in /atoms/
c  ity  - number of atom types among neighbours
c    do i=1,ity  ! for each neighbour type
c      natnum(i) - atomic number of this type 
c      ndt(i)   - number of neighbours of this type
c    enddo 
c  iring = 2 - atom belongs to 2 planar rings
c        = 1 - atom belong to 1 planar ring
c        = 0 - atom does not belong to a ring or ring is not planar or
c              it belongs to several rings
c  asym(ik) - symmetry label (converted later to upper case)
c
c for H-atoms ONLY
c      ihp(ik)   - atom number in periodic table of parent atom
c      isp(ik)   - number of cov-bonded nei to parent atom
c      ihnei(ik) - number of cov-bonded H-atom to parent atom
c
c.....exception - group type
      character group(nat)*40
      common /except/ iexp(nat), group
c     local params 
      parameter(maxneityp=10) ! maximum neighbour types
      dimension jatum(maxneityp), jndt(maxneityp)
c      dimension matum(maxneityp), mndt(maxneityp)
      dimension rdum(100)
      character acom*80, line*200, abank*60
      common /mpoles1/ ifind, Pv, sPv, rkap(2), Plm(0:4,-4:4)
      common /abank/ abank
      logical lPG, lPGi
      common /lPG/ lPG
      parameter( maxchar = 40 ) ! max number of character variables in the line
      character vchar(maxchar)*8
      integer istar(maxchar)
      common /line_extract/ nch,istar,vchar
      common /iprint/ iprint    
      common /ldb/ ldb, ntot, nfound, nmissed
      data zero /0.00D00/
c     atomic number
      ityp = itype(ik)
c     redo the search for neighbors        
      iprint = 0
      call FindNei( ik, zero )
c     initialize arrays
      line   = ' '
      acom   = ' '
      groupj = ' '
      sPv = 0.01D00
      Pv = Pval(ityp)
      asymr = ' '
      asymi = ' '
      rkap(1)=1.00D00
      rkap(2)=1.00D00
      do l=0,4
        do m=-l,l
          Plm(l,m)=zero
        enddo
      enddo
c
      asymi = asym(ik)
      call UpCase(asymi)
c
      write(io,'(/a)') ' Searching DataBank for a match : '
      write(io,*) ' ityp  = ',ityp
      write(io,*) ' iring = ',iring
      write(io,*) ' inei  = ',inei
      write(io,*) ' ity   = ',ity
      write(io,*) ' imemb = ',imemb
      write(io,'(3a)')'  asymi = [',asymi,']'
      write(io,*) ' lPG   = ',lPG  
      do i=1,ity
        write(io,*)' type = ',i,' atom = ',natnum(i),' neib = ',ndt(i)
      enddo
      write(io,'(a)') ' '
c     
      ifind = 0      
c     rewind databank
      rewind(idb)
c
c     reading DataBank
c
1     read(idb,99,end=100) line
c       find atom record
        if(line(1:4).ne.'ATOM') goto 1
c        write(io,99) ' Reading ATOM = ',line
c       find atomic number
        call ReadChar(line,2,i1,i2)
        read(line(i1-1:i2+1),*) janum
c        write(io,*) ' atom number  = ',janum
        if(janum.ne.ityp) goto 1
        if(iverbose.ge.2) write(io,99) ' Reading atom = ',line
c       save comment
        acom=' '
        call ReadChar(line,4,i1,i2)
        ijunk=ils(line)-i1+1
        acom(1:ijunk)=line(i1:ils(line))
c       assign defaults for lPGi (same as for lPG in DefLoc)   
        lPGi = .true.
        if( inei.ge.4 ) lPGi = .false.
c
c.......H-atoms - special case
c
        if(ityp.eq.1) then
110       call ReadFile(idb,line)
          if(vchar(1).ne.'PRNT') goto 110
          read(line,*) jhp,jsp,jhnei
          asymr = vchar(5)
c         check for EXCEPTIONS!!!!
          call ExceptDB(ik,groupj,*1)
          if(iverbose.ge.2) write(io,112) 
     *    'jhp =',jhp,'jsp =',jsp,'jhnei=',jhnei,'asymr=[',asymr,']'
          if( jhp  .ne.ihp(ik)   ) goto 1 
          if( jsp  .ne.isp(ik)   ) goto 1 
          if( jhnei.ne.ihnei(ik) ) goto 1 
          if( asymr.ne.asymi     ) goto 1  
          ifind=1
          goto 12
        endif
112     format(2x,3(1x,a,i3),1x,3a,1x,a,l1)
c
c.......all other atoms - general case
c
c       find NEIG line
11      call ReadFile(idb,line)
        if(vchar(1).ne.'NEIG') goto 11
        asymr = vchar(6)
        read(line,*) jnei, jty, jring, jmemb
        if(vchar(8).eq.'F') lPGi = .false.
        if(vchar(8).eq.'T') lPGi = .true.
        if(iverbose.ge.2) 
     *    write(io,112)  'jnei=',jnei, 'jty=',jty, 'jring=',jring, 
     *                   'asymr=[',asymr,']', 'planar=',lPGi,
     *                   'jmemb=',jmemb
c       skip if number of neighbors or types not the same, 
c         or iring does not match, or symmetry is different
        if(  inei.ne.jnei  ) goto 1
        if(   ity.ne.jty   ) goto 1
        if( iring.ne.jring ) goto 1
        if( imemb.ne.jmemb ) goto 1
c        if( imemb.ne.jmemb ) goto 1  
c       if there is no symmetry label for this atom in databank,
c       i.e. asymr = ' ', then do not make the check for symmetry
c       consistency because we need to do it only for special cases anyway
        if( asymr.ne.' ' ) then 
           if( asymr.ne.asymi ) goto 1
        endif   
c
        if( lPG.neqv.lPGi ) goto 1
c       now check the neighbour types 
c       (actual types and number of nei of each type)
        do j=1,jty
          call ReadFile(idb,line)
          read(line,*) ijunk, jatum(j), jndt(j)
          if(iverbose.ge.2) write(io,4889) j, jatum(j), jndt(j)
4889      format(3x,'type =',i3,' jatum =',i3,' jndt =',i3)
        enddo                
        im=0
        do i=1,ity
         im=0
         do j=1,jty
           if(natnum(i).eq.jatum(j)) then
              if(ndt(i).eq.jndt(j)) im=1
           endif
         enddo 
         if(im.eq.0) goto 1 
        enddo
c       check for EXCEPTIONS!!!!
        call ExceptDB(ik,groupj,*1)
c
c       ok, let's think now :
c       atomic number - passed
c       number of neighbours - passed
c       number of neighbour types - passed        
c       each of neighbour types - passed
c       number of nei of each type - passed
c       local symmetry - passed
c       if exception was requested - passed    
c       well, looks like we found our atom
c
        ifind=1
c       read Pv & kappas
12      call ReadFile(idb,line)
        if(vchar(1).ne.'PVAL') goto 12
        read(line,*) Pv, rkap, sPv 
c       read multipoles until reach another 'ATOM' record
13      call ReadFile(idb,line)
          if(vchar(1).eq.'ATOM') goto 100
          if(vchar(1).ne.'PLMS') goto 100
          call DNumLine(line,ir)
          read(line,*) (rdum(k),k=1,ir)
          ij=1
          do i=1,ir/3
            l=nint(rdum(ij))
            ij=ij+1
            m=nint(rdum(ij))
            ij=ij+1
            plm(l,m)=rdum(ij)
            ij=ij+1
          enddo 
        goto 13
c
c     end loop reading DataBank
c      
100   continue        
c     output
      n = ils(atom(ik))
      ntot = ntot + 1
      if(ifind.eq.0) then
        write(io,60) atom(ik)(:n), ' NOT '
        nmissed = nmissed + 1
      else if(ifind.eq.1) then
        write(io,60) atom(ik)(:n), ' '
        write(io,62) atom(ik)(:n), acom(1:ils(acom))
        if( iexp(ik).ne.0 ) write(io,63) groupj(:ils(groupj))
        nfound = nfound + 1
      else
        call Err('ReadDB','Unknown value of ifind')  
      endif  
      write(io,61) atom(ik)(:n), Pv,rKap, sPv, ((Plm(l,m),m=-l,l),l=0,4)
60    format(/' ReadDB :: atom ',a,' was',a,'found in the DataBank')
62    format(/' DB comment for ',a,' = [',a,']')
63    format(/' Functional group = [',a,']')
61    format(/' Pseudoatom parameters assigned to ',a,' : ',//,
     * 3x,'Pval, Kappa, Kappa'', sigmaPv = ',4f12.4,/,
     * 3x,' 2nd Monopole:',f7.3,/,
     * 3x,'      Dipoles:',3f7.3,/,
     * 3x,'  Quadrupoles:',5f7.3,/,
     * 3x,'    Octupoles:',7f7.3,/
     * 3x,'Hexadecapoles:',9f7.3) 
c
99    format(100a)
      return
c
101   call Err('ReadDB','Error reading DataBank file...')
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c reads line from file inp, skips commented lines (!)
c converts line to upper case and 
c calls read_line subroutine
      subroutine ReadFile(inp,line)
      character*(*) line, ajunk*80
1     read(inp,'(a)',end=2) line
        if(line(1:1).eq.'!') goto 1
      call UpCase(line)
      call ReadLine(line)
      return
2     ajunk = ' '
      write(ajunk,*) 'Error reading unit = ',inp
      call Err('ReadFile',ajunk)
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c handles exception treatment when reading db
      subroutine ExceptDB(ik,groupj,*)
      include 'lsdb.inc'
      character line*200, groupj*(*)
      parameter( maxchar = 40 ) ! max number of character variables in the line
      character vchar(maxchar)*8
      integer istar(maxchar)
      common /line_extract/ nch,istar,vchar
      character group(nat)*40
      common /except/ iexp(nat), group
      common /files/ imas, io, ixp, idb, icon
      common /verbosity/ iverbose
c
601   format(/' Functional group = [',a,'], but no functional group',
     *      ' specified for this entry..',/,' Skipping DB entry...'/)  
602   format(3x,'Exception entry in DataBank :')
603   format(' Exception was not identified for this atom - ',
     *       'Skip DB entry...',/)
604   format(6x,'functional group   = [',a,']')
605   format(6x,'group read from DB = [',a,']')
c
      call ReadFile(idb,line)
      if( vchar(1).ne.'EXCEPT' ) then
        backspace(idb)
c       so there is no exception flag in db, but now we have to check 
c       whether the current atom is an exception
        if( iexp(ik).ne.0 ) then
          if(iverbose.ge.2) write(io,601) group(ik)(:ils(group(ik)))
          return 1
        endif          
c       ok, no exception for atom and no exception in db - go on...
        return
      endif   
c     ok, we have an exception - read functional group
      if(iverbose.ge.2) write(io,602)
c     if exception was not identified to this atom, 
c     why should we continue reading the databank ?
      if( iexp(ik).eq.0 ) then
        if(iverbose.ge.2) write(io,603)
        return 1
      endif  
c     ok, continue
      backspace(idb)
      read(idb,'(a)') line
      call upcase(line)
      i1 = index(line,'GROUP')
      i2 = ifs( line( i1:len(line) ) )
      i3 = i1 + 5 + i2 
      groupj = line( i3:ils(line) )
      if( iverbose.ge.2 ) then
        write(io,604) group(ik)(:ils(group(ik)))
        write(io,605) groupj(:ils(groupj))
      endif  
      if( groupj.ne.group(ik) ) then
        if( iverbose.ge.2 ) 
     *    write(io,'(a/)') '   Wrong functional group... :( '        
        return 1
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c rescales Pv populations of fragment(s) if requested
c
      subroutine rescalePv
      include 'lsdb.inc'
      character line*255, ans*1
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym    
      common /atoms1/ iut(nat), isft(nat), lmx(nat), occ(nat), iaf(nat) 
      common /mpoles/ Pv(nat), sPv(nat), Plm(nat,0:4,-4:4), iar(nat,11)
      logical lring
      character fname(2)*60
      common /args/ fname, lring, isca, sumofPV
      common /files/ imas, nout, ixp, idb, icon
      data zero /0.00D00/
      dimension iflag(nat), idum(nat), chfrag(nat)
      logical lall, lcfrag
      parameter( maxchar = 40 ) ! max number of character variables in the line
      character vchar(maxchar)*8
      integer istar(maxchar), ise(maxchar,2)
      common /line_extract/ nch,istar,vchar
c
      line = ' '
c
c     method of rescaling
c
c     if isca was already specified in input file....
      if( isca.ne.0) goto 333
      isca = 3
      ans=' '
      write(*,97) 
     *  ' Use difference (d), fraction (f) or sigmas (s) [s] : '
      read(*,99) ans
      if(ans.eq.'d'.or.ans.eq.'D') isca=1         
      if(ans.eq.'f'.or.ans.eq.'F') isca=2        
333   if(isca.eq.1) write(nout,430) 'difference'
      if(isca.eq.2) write(nout,430) 'fraction'
      if(isca.eq.3) write(nout,430) 'sigmaPv'
430   format(/' Using ',a,' method for rescaling...')           
c
      nfrag = 0
      lall = .true.
      lcfrag = .false.
      call rzero(chfrag,nat)
c
c     first we try to read fragment info from lsdb.inp file
c     (same syntax as for INTEREN in xd.mas)
c
      ij=45
      open(ij,file='lsdb.inp',form='formatted',status='old',err=202)
200   call ReadInputLine(ij,line,*201)     
        call upcase(line)
        if(line(1:5).eq.'NFRAG') read(line(6:),*) nfrag
        if(line(1:4).eq.'FRAG') lall = .false.
        if(line(1:6).eq.'CHFRAG') then
           lcfrag = .true.
           call ReadLine(line)
           call DNumLine(line,ndum)
           read(line,*) (chfrag(i),i=1,ndum)
        endif   
      goto 200       
201   continue
      if(nfrag.eq.0) goto 202      
      if(lall) goto 202      
c.....ok, we just read number of fragments from lsdb.inp file
      goto 1002
c
c     either lsdb.inp was not found, or fragment info was not found in lsdb.inp file
c     use all atoms in fragment 1
c
202   continue
      nfrag = 1
      write(   *,203)
      write(nout,203)
203   format(/' Using ALL atoms for fragment 1 !')
c
c.....loop over fragments
c
1002  continue
      do 1000 k=1,nfrag  
c.......open filename for k-th fragment
        write(   *,3889) k
        write(nout,3889) k
3889    format(/' ***** Rescaling Pv''s in fragment = ',i3,' *****')
c.......if all atoms belong to 1st fragments
        if( lall ) then
          do i=1,nat
            iflag(i) = 1            
          enddo        
          goto 47
        endif
c.......weird way of reading atom numbers of k-th fragment
        call izero(iflag,nat)
        line=' '
        rewind(ij)
46      call ReadInputLine(ij,line,*47)
        call upcase(line)
        if(line.eq.' ') goto 46
        if(line(1:4).ne.'FRAG') goto 46
        line(1:4)='    '
        call DNumLine(line,inum)  
c        write(*,*) ' k, inum = ',k,inum      
        read(line,*,err=1005) (idum(j),j=1,inum)
        if( idum(1).ne.k ) goto 46
        do i=2,inum
          j = idum(i)
          if(j.gt.0) iflag(j) = 1
          if(j.lt.0) then
            istart=idum(i-1)
            iend=abs(j)
            do m=istart,iend
              iflag(m) = 1
            enddo
          endif
        enddo
        goto 46                     
47      continue
c.......pre-calculate       
        Ptot = zero
        Pfre = zero
        Psig = zero
        do 100 i=1,ia
          if( iflag(i).eq.0 ) goto 100
          Ptot = Ptot + Pv(i)
          Pfre = Pfre + Pval( itype(i) )*occ(i)
          Psig = Psig + sPv(i)
100     continue
        delP = Pfre - Ptot
        dumP = abs(delP)/Pfre*1.D02
        write(   *,3890) Pfre, Ptot, delP, dumP
        write(nout,3890) Pfre, Ptot, delP, dumP
3890    format(/' Sum-of-Pv[free]     = ',f12.6,' el.',/,
     *          ' Sum-of-Pv[databank] = ',f12.6,' el.',/,
     *          ' Difference          = ',f12.6,' el. ( ',f6.2,' % )')     
        if( isca.eq.3 ) then    
          write(   *,3891) Psig
          write(nout,3891) Psig
3891      format(' Sum of sigma(Pv)    = ',f12.6)
        endif
c.......new number of valence electrons
c       charge of fragment was NOT read from lsdb.inp file
        if( .not. lcfrag ) then
          write(*,3892) k
3892      format(' Enter net charge of fragment ',i3,' [0.] : ',$)
          read(*,99) line
          if(line.ne.' ') read(line,*) chfrag(k)
        endif   
        Pres = Pfre - chfrag(k) 
        write(   *,3893) Pres, chfrag(k)  
        write(nout,3893) Pres, chfrag(k)
3893    format(/' Rescaling Sum-of-Pv to ',f13.6,' el.',3x,
     *          '( Charge = ',sp,f8.3,' el. )')
c.......perform actual rescaling          
        Pnew = zero
        if(isca.eq.1) write(nout,3894) '+'
        if(isca.eq.2) write(nout,3894) '*'
        if(isca.eq.3) write(nout,3894) 's'
3894    format(/4x,'NO  ATOM',7x,'Occup.    Pv(DB)      Pv(free)',
     *             '    ',a,'Corr       Pv(new)'/)
c.......difference (isca=1)
        if(isca.eq.1) then
          delta = (Pres-Ptot) / ia
          do 101 j=1,ia
            if( iflag(j).eq.0 ) goto 101
            Pv1 = Pv(j) + delta
            write(nout,387) j, atom(j), occ(j), Pv(j), 
     *                      Pval(itype(j))*occ(j), delta, Pv1
            Pv(j) = Pv1
            Pnew = Pnew + Pv(j)
101       continue
          write(nout,385) Ptot,Pfre,Pres-Ptot,Pnew
c.......fraction (isca=2)
        else if(isca.eq.2) then 
          delta = Pres / Ptot
          do 102 j=1,ia
            if( iflag(j).eq.0 ) goto 102
            Pv1 = Pv(j)*delta
            write(nout,387) j, atom(j), occ(j), Pv(j), 
     *                      Pval(itype(j))*occ(j), delta, Pv1
            Pv(j) = Pv1
            Pnew = Pnew + Pv(j)
102       continue
          write(nout,385) Ptot,Pfre,Pres-Ptot,Pnew
c.......based on sigmas (isca=3)
        else if(isca.eq.3) then 
          do 103 j=1,ia
            if( iflag(j).eq.0 ) goto 103
            delta = ( (Pres-Ptot) / Psig ) * sPv(j)
            Pv1 = Pv(j) + delta
            write(nout,387) j, atom(j), occ(j), Pv(j), 
     *                      Pval(itype(j))*occ(j), delta, Pv1
            Pv(j) = Pv1
            Pnew = Pnew + Pv(j)
103       continue
          write(nout,385) Ptot,Pfre,Pres-Ptot,Pnew
        endif
c
c.....end loop over fragments
c
1000  continue
      close(ij)
      return
c
1005  call Err('rescalePv','Error reading fragments from lsdb.inp')
c
385   format(/20x,'SUM ',4f12.5)
387   format(i6,2x,a,1x,f7.3,2f12.5,sp,f12.5,ss,f12.5)
97    format( a,$)
98    format(/a,$)
99    format(100a)  
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c transforms CARTESIAN<-->FRACTIONAL 
c iflag=0 from fractional {x1,y1,z1} to cartesian  {x2,y2,z2}
c iflag=1 from cartesian  {x1,y1,z1} to fractional {x2,y2,z2} 
      subroutine rabs(iflag,x1,y1,z1,x2,y2,z2)
      common /transf/ rm(3,3), rm1(3,3)
      if(iflag.ne.0.and.iflag.ne.1) call Err('rabs','iflag is wrong...')
c     fractional-->cartesian
      if(iflag.eq.0) then
        x2 =  rm(1,1)*x1 +  rm(1,2)*y1 +  rm(1,3)*z1
        y2 =  rm(2,1)*x1 +  rm(2,2)*y1 +  rm(2,3)*z1
        z2 =  rm(3,1)*x1 +  rm(3,2)*y1 +  rm(3,3)*z1
c     cartesian-->fractional
      else
        x2 = rm1(1,1)*x1 + rm1(1,2)*y1 + rm1(1,3)*z1
        y2 = rm1(2,1)*x1 + rm1(2,2)*y1 + rm1(2,3)*z1
        z2 = rm1(3,1)*x1 + rm1(3,2)*y1 + rm1(3,3)*z1
      endif 
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ReadLine  taken from HIPPO98 (Z. Su and A. Volkov). 
c This subroutine reads a line of arbitrary length. It finds our whether n-th 
c character parameter in the line is marked with (*). If yes, then istar(n)=1, 
c if not then istar(n)=0. The keywords are stored in "vchar" array. 
c Then this subroutine replaces the character variables with blanks so the numerical
c variables can be read easily in free format
c "nch" is the number of character variables in the line
c
      Subroutine ReadLine(line)
c
      parameter( maxchar = 40 ) ! max number of character variables in the line
      character line*(*), vchar(maxchar)*8, aform*10
      integer istar(maxchar), ise(maxchar,2)
      common /line_extract/ nch,istar,vchar
      
cAV 11/16/00  ignore first 4 characters if it's a keyword      
ccccc      if(line(1:1).ne.' ')  line(1:4)='    '
cendAV
      Do i=1,maxchar
         istar(i) = 0
         vchar(i) = ' '
      End Do

      nch  = 0
      
      call DNumLine(line,inum)
c      write(*,*) ' inum = ',inum
      do i=1,inum
        call ReadChar(line,i,i1,i2)
        aform = ' '
        write(aform(1:),38) i2-i1+1
38      format('(f',i4,'.0)')
        call RmBlanks(aform)
c        write(*,'(100a)') ' line(i1:i2) = [',line(i1:i2),']    ',
c     *                    ' aform = ',aform
c.......decide if i-th variable is real/integer or character
        read(line(i1:i2),aform,err=104) rdum
c        write(*,*) ' rdum = ',rdum
        goto 105
c.......character
104     nch = nch + 1
c        write(*,*) '   nch = ',nch
        if(line(i1:i1).eq.'*') then
          istar(nch)=1
          vchar(nch)=line(i1+1:i2)  
          call upcase(vchar(nch)) 
        else
          vchar(nch)=line(i1:i2)
          call upcase(vchar(nch))
        endif 
        ise(nch,1)=i1
        ise(nch,2)=i2       
105     continue
      enddo 
c.....now clean the line of character variables
      do i=1,nch
        line(ise(i,1):ise(i,2))=' '
      enddo
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c returns serial numbers of 1st (i1) and last(i2) symbols of the 
c is-th character (CP) parameter in character variable "line" of arbitrary length
c 11/14/00: The length of the "line" is determined with built-in function "len"
c 07/24/01: Fixed bug when the first character parameter starts at the 1st column
c 06/26/01: Almost completely rewritten 
      subroutine ReadChar(line,is,i1,i2)
      character*(*) line
c
      i1=0
      i2=0
      ns=len(line)
c
c49      write(49,'(/a/3a)') ' New call to READCHAR Line : ','[',
c49     * line(:ils(line)),']'
      js=0
      do i=1,ns-1
        if(line(i:i).ne.' '.and.line(i+1:i+1).eq.' ') js=js+1  
c49        write(49,567) i,i,line(i:i),i+1,i+1,line(i+1:i+1),js     
c49567     format(4x,2(' line(',i3,':',i3,') = ',a,3x),'js=',i2)
        if(js.eq.is) then
          i2=i
c49          write(49,*) ' READCHAR :: is,i2 = ',is,i2
          do j=i2,1,-1
             if(line(j:j).eq.' '.and.line(j+1:j+1).ne.' ') then 
               i1=j+1
               return
             endif
          enddo 
          i1=1
          return
        endif
      enddo
      call Err('ReadChar','Unexpected STOP in ReadChar...')
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c number of variables in line, separated by al least 1 space
      subroutine DNumLine(line,inum)
      character line*(*)
      llen=len(line)
      inum=0
      do i=1,llen-1
        if(line(i:i).ne.' '.and.line(i+1:i+1).eq.' ') inum=inum+1
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c removes blanks from a character variable
      subroutine RmBlanks(op)
      Character*(*) op
      Character*256 line
      line= ' '
      if( len(op) .gt. len(line) ) 
     *   call Err('RmBlanks','len(op).gt.len(line)')
      ij=0  
      do i=1,len(op)
        if(op(i:i).ne.' ') then
          ij=ij+1
          line(ij:ij)=op(i:i)
        endif
      enddo      
      op=' '
      op(1:len(op))=line(1:len(op))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c convert string to lowercase
      Subroutine LoCase(string)
      Character*(*) string
      Do 1000 i=1,Len(string)
         If ( LLE(string(i:i),'Z') .And. LGE(string(i:i),'A') )
     +        string(i:i) = Char( IChar(string(i:i)) + 32 )
 1000 Continue
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c convert string to uppercase
      Subroutine UpCase(string)
      Character*(*) string
      Do 1000 i=1,Len(string)
         If ( LLE(string(i:i),'z') .And. LGE(string(i:i),'a') )
     +        string(i:i) = Char( IChar(string(i:i)) - 32 )
 1000 Continue
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c error handling
      subroutine Err(namz,mzss)
      character*(*) namz
      character*(*) mzss
      common /files/ imas, nout, ixp, idb, icon
      im = ils(mzss)
      write(nout,1) namz, mzss(:im)
      write(   *,1) namz, mzss(:im)
1     format(/' ERROR *** ',a,' *** ',a,/)
      close(nout)
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rzero(r,n)
      dimension r(n)
      data zero /0.00D00/
      do i=1,n
        r(i)=zero
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine izero(i,n)
      dimension i(n)
      do j=1,n
        i(j)=0
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            BLOCK DATA
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     covalent radii, atomic symbols, rm & rm1
      block data
c
      common /transf/ rm(3,3), rm1(3,3)
      data rm  /1.,0.,0.,0.,1.,0.,0.,0.,1./
      data rm1 /1.,0.,0.,0.,1.,0.,0.,0.,1./
c
      common /ifra/ ifra
      data ifra /0/
c
      common /covrad/ covrad(54)
      data covrad/      
     *0.40,0.93,
     *1.23,0.90,0.82,0.77,0.75,0.73,0.72,0.71,
     *1.54,1.36,1.18,1.11,1.06,1.02,0.99,0.98,
     *2.03,1.74,1.44,1.32,1.22,1.18,1.17,1.17,1.16,1.15,1.17,1.25,1.26,
     *                                        1.22,1.20,1.16,1.14,1.89,
     *2.16,1.91,1.62,1.45,1.34,1.30,1.27,1.25,1.25,1.28,1.34,1.41,1.44,
     *                                        1.41,1.40,1.36,1.33,1.31/
c
      character symbat(92)*2
      common /atsym/ symbat
      data symbat/        'H ','He','Li','Be','B ','C ','N ',
     1'O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar',
     2'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     3'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     4'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',
     5'Tr','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     6'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
     7'W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
     8'At','Rn','Fr','Ra','Ac','Th','Pa','U '/
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            FUNCTIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     finds the length of the vector x
      function RLength(x)
      dimension x(3)
      data zero /0.00D00/
      rlength = zero
      do i=1,3
        rlength=rlength+x(i)**2
      enddo
      rlength=sqrt(rlength)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c returns the number of atom to which we should go from il-th atom if 
c we came from ik-th atom
c if ichoice = 0 - nowhere to go....
      function ichoice(ik,il)
      include 'lsdb.inc'
      common /conmat/ icmat(nat,nbond), nba(nat)
      common /choices/ nbaf(nat,nbond)   
      common /files/ imas, nout, ixp, idb, icon
      dimension ird(0:maxmem)
      ichoice = 0
c      write(nout,'(/a,2i4)') 'ik --> il = ',ik,il
c     loop over neigbors of ik-th atom
      do i=1,nba(il)
         ipath = nbaf(il,i)
         iatf = icmat(il,i)
         if( ipath.ne.0 .and. iatf.ne.ik ) goto 51
c         if( ipath.ne.0 ) goto 51
      enddo
      return
c
51    continue
      ichoice = iatf
c     "close" this path (1->2)
      nbaf(il,i) = 0
c      write(nout,*) ' closing path il <-> i : ', il, ' <-> ', iatf
c     also close this path from the other "side" (2->1)
c      do j=1,nba(iatf)
c         if( icmat(iatf,j) .eq. il ) nbaf(iatf,j) = 0
c      enddo 
      return
      end    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c returns 1 if found ring is already in the list
c         0 if this is a new ring
      integer function ichk_ring(ir,id)   
      include 'lsdb.inc'
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      dimension id(maxmem)
      ichk_ring=0
      do 900 i=1,nor
        if( ir.ne.noar(i) ) goto 900
        im = 0
        do j=1,noar(i)
          do k=1,ir
            if(nra(i,j).eq.id(k)) im = im + 1
          enddo 
        enddo
        if( im.eq.noar(i) ) then
          ichk_ring=1
          return
        endif  
900   continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c very crude: number of valence electrons based on atomic number
      function Pval(iat)
      data zero /0.00D00/
      Pval = zero
      if(iat.le. 2)               Pval = real(iat)
      if(iat.ge. 3.and.iat.le.10) Pval = real(iat- 2)
      if(iat.ge.11.and.iat.le.18) Pval = real(iat-10)
      if(iat.ge.19.and.iat.le.36) Pval = real(iat-18)
      if(iat.ge.37.and.iat.le.54) Pval = real(iat-36)
      if(iat.ge.55.and.iat.le.86) Pval = real(iat-54)
      if(iat.ge.87)               Pval = real(iat-86)
      return      
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cKNJ
c compares two NON-HYDROGEN atoms for chemical equivalency 
c if iSimple <> 0 : compare just the number of neighbors. that's all...
c if iCs <> 0 : also compares local symmetries
c if iCg <> 0 : also compares functional groups 
      logical function SameAtoms( ia1, ia2, iCs, iCg, iSimple, iRg )
      include 'lsdb.inc'
      logical pla, pla2, ir1, plan, plan2, ir2, cyc1, cyc2
      character line*40
      common /files/ imas, nout, ixp, idb, icon 
      logical planar
      common /rings/ nor, noar(maxrings), nra(maxrings,maxmem), 
     *               planar(maxrings), ipldum(maxrings)
      common /iprint/ iprint  
      character atom(nat)*8, asym(nat)*5
      common /atoms/ ia, atom, xf(nat,3), xc(nat,3), itype(nat), asym 
      character group(nat)*40
      common /except/ iexp(nat), group
c     nnei - number of neighbors
c     jnei - sequence # of each neighbor in  /atoms/ 
c     rnei - distance to each neigbor
c     jtnei - number of atoms types among neighbors 
c     janei - atomic number of each neighbor type
c     jnnei - number of neighbors of each type
c     jsnei - sequence of neighbors in /atoms/ for each type 
      common /neighbors/ nnei(nat), jnei(nat,nbond), rnei(nat,nbond),
     *                   jtnei(nat), janei(nat,nbond), jnnei(nat,nbond),
     *                   jsnei(nat,nbond,nbond)        
      common /verbosity/ iverbose
c
      SameAtoms = .false.
c
      ma1 = 0
      ma2 = 0
      i1 = ils( atom(ia1) )
      i2 = ils( atom(ia2) )
      if(iverbose.ge.3) write(nout,528) atom(ia1)(:i1), atom(ia2)(:i2)
528   format(/' Checking if atoms ',a,' and ',a,
     *        ' are chemically equivalent:')     
c.....
529   format(3x,'atoms are NOT chemically equivalent : ',a)
c.....different atoms at all
      if( itype(ia1) .ne.itype(ia2) ) then
        if(iverbose.ge.3) write(nout,529) ' Diff. atomic numbers'
        return
      endif  
c.....test the equivalency
      if( nnei(ia1).ne.nnei(ia2) ) then
        if(iverbose.ge.3) write(nout,529) ' Diff. number of neighbors'
        return ! different number of neighbors
      endif   
c.....for simple comparison -> stop here      
      if( iSimple.ne.0 ) goto 10
c.....      
      if( jtnei(ia1).ne.jtnei(ia2) ) then
        if(iverbose.ge.3) write(nout,529) 
     *   ' Diff. number of neighbor types'
        return ! different number of neighbor types
      endif  
      do i=1,jtnei(ia1)
        ityp = janei(ia1,i)
        jtyp = janei(ia2,i)
        if( ityp.ne.jtyp ) then
          line = ' '
          write(line,'(a,i2,a)')
     *     ' Diff. atomic numbers of ',i,'-th type'
          if(iverbose.ge.3) write(nout,529) line
          return ! different atomic number of i-th neighbor type
        endif
        nai = jnnei(ia1,i)
        naj = jnnei(ia2,i)
        if( nai.ne.naj ) then
          line = ' '
          write(line,'(a,i2,a)')
     *     ' Diff. num. of nei. of ',i,'-th type'
          if(iverbose.ge.3) write(nout,529) line
          return ! different number of neighbors of i-th neighbor type
        endif
      enddo
c.....check number of ring members
      if( iRg.ne.0 ) then
      call IsItCyclic(ia1, cyc1, pla, pla2, ir1, m1)
           ma1= ma1+m1
      call IsItCyclic(ia2, cyc2, plan, plan2, ir2, m2)
           ma2= ma2+m2
      if( iRg.ne.0 .and. ma1.ne.ma2 ) then
        if(iverbose.ge.3) write(nout,529) ' Diff. ring memb numbers'
        return  ! different local symmetry
      endif
      endif
c
10    continue
c.....check the local symmetry
      if( iCs.ne.0 .and. asym(ia1).ne.asym(ia2) ) then
        if(iverbose.ge.3) write(nout,529) ' Diff. local symmetry'
        return  ! different local symmetry
      endif
c.....check functional group
      if( iCg.ne.0 .and. group(ia1).ne.group(ia2) ) then
        if(iverbose.ge.3) write(nout,529) ' Diff. functional groups'
        return  ! different local symmetry
      endif
c
      SameAtoms = .true.
      if(iverbose.ge.3) 
     * write(nout,*) '  YES, these atoms are chemically equivalent !'
      return      
      end
cKNJ
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function asind(arg)
      parameter( pi=3.1415926535897932384626433832795028841971694D0 )
      parameter( todeg = 180.00D00/pi )
      parameter( torad = pi/180.00D00 )
      asind = todeg*asin(arg)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function sind(arg)
      parameter( pi=3.1415926535897932384626433832795028841971694D0 )
      parameter( todeg = 180.00D00/pi )
      parameter( torad = pi/180.00D00 )
      sind = sin(arg*torad)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function cosd(arg)
      parameter( pi=3.1415926535897932384626433832795028841971694D0 )
      parameter( todeg = 180.00D00/pi )
      parameter( torad = pi/180.00D00 )
      cosd = cos(arg*torad)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c distance between two atoms in cartesian coodinates
      function Dist(x1,y1,z1,x2,y2,z2)
      data zero /0.00D00/
      Dist = zero
      dx = x2 - x1
      dy = y2 - y1
      dz = z2 - z1
      Dist = dx*dx + dy*dy + dz*dz
      Dist = sqrt( Dist )     
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c returns atomic number or 0 if not recognized
      function iAtNum(atom)
      character*(*) atom
      character at*2
      at=' '
      at(1:2)=atom(1:2)
      call UpCase(at)
      iAtNum=0 ! atomic number
      if(at(1:1).eq.'H')  iAtNum =   1
      if(at(1:1).eq.'B')  iAtNum =   5      
      if(at(1:1).eq.'C')  iAtNum =   6
      if(at(1:1).eq.'N')  iAtNum =   7
      if(at(1:1).eq.'O')  iAtNum =   8
      if(at(1:1).eq.'F')  iAtNum =   9
      if(at(1:1).eq.'P')  iAtNum =  15
      if(at(1:1).eq.'S')  iAtNum =  16
      if(at(1:1).eq.'K')  iAtNum =  19
      if(at(1:1).eq.'V')  iAtNum =  23
      if(at(1:1).eq.'Y')  iAtNum =  39
      if(at(1:1).eq.'I')  iAtNum =  53
      if(at(1:1).eq.'W')  iAtNum =  74
      if(at(1:1).eq.'U')  iAtNum =  92
c
      if(at(1:2).eq.'HE') iAtNum =   2
      if(at(1:2).eq.'LI') iAtNum =   3 
      if(at(1:2).eq.'BE') iAtNum =   4 
      if(at(1:2).eq.'NE') iAtNum =  10 
      if(at(1:2).eq.'NA') iAtNum =  11 
      if(at(1:2).eq.'MG') iAtNum =  12 
      if(at(1:2).eq.'AL') iAtNum =  13 
      if(at(1:2).eq.'SI') iAtNum =  14 
      if(at(1:2).eq.'CL') iAtNum =  17 
      if(at(1:2).eq.'AR') iAtNum =  18 
      if(at(1:2).eq.'CA') iAtNum =  20 
      if(at(1:2).eq.'SC') iAtNum =  21 
      if(at(1:2).eq.'TI') iAtNum =  22 
      if(at(1:2).eq.'CR') iAtNum =  24 
      if(at(1:2).eq.'MN') iAtNum =  25 
      if(at(1:2).eq.'FE') iAtNum =  26 
      if(at(1:2).eq.'CO') iAtNum =  27 
      if(at(1:2).eq.'NI') iAtNum =  28 
      if(at(1:2).eq.'CU') iAtNum =  29  
      if(at(1:2).eq.'ZN') iAtNum =  30 
      if(at(1:2).eq.'GA') iAtNum =  31 
      if(at(1:2).eq.'GE') iAtNum =  32 
      if(at(1:2).eq.'AS') iAtNum =  33 
      if(at(1:2).eq.'SE') iAtNum =  34 
      if(at(1:2).eq.'BR') iAtNum =  35 
      if(at(1:2).eq.'KR') iAtNum =  36 
      if(at(1:2).eq.'RB') iAtNum =  37 
      if(at(1:2).eq.'SR') iAtNum =  38 
      if(at(1:2).eq.'ZR') iAtNum =  40 
      if(at(1:2).eq.'NB') iAtNum =  41 
      if(at(1:2).eq.'MO') iAtNum =  42 
      if(at(1:2).eq.'TC') iAtNum =  43 
      if(at(1:2).eq.'RU') iAtNum =  44 
      if(at(1:2).eq.'RH') iAtNum =  45 
      if(at(1:2).eq.'PD') iAtNum =  46 
      if(at(1:2).eq.'AG') iAtNum =  47 
      if(at(1:2).eq.'CD') iAtNum =  48 
      if(at(1:2).eq.'IN') iAtNum =  49 
      if(at(1:2).eq.'SN') iAtNum =  50 
      if(at(1:2).eq.'SB') iAtNum =  51 
      if(at(1:2).eq.'TE') iAtNum =  52 
      if(at(1:2).eq.'XE') iAtNum =  54 
      if(at(1:2).eq.'CS') iAtNum =  55 
      if(at(1:2).eq.'BA') iAtNum =  56 
      if(at(1:2).eq.'LA') iAtNum =  57 
      if(at(1:2).eq.'CE') iAtNum =  58 
      if(at(1:2).eq.'PR') iAtNum =  59 
      if(at(1:2).eq.'ND') iAtNum =  60 
      if(at(1:2).eq.'PM') iAtNum =  61 
      if(at(1:2).eq.'SM') iAtNum =  62 
      if(at(1:2).eq.'EU') iAtNum =  63 
      if(at(1:2).eq.'GD') iAtNum =  64 
      if(at(1:2).eq.'TB') iAtNum =  65 
      if(at(1:2).eq.'DY') iAtNum =  66 
      if(at(1:2).eq.'HO') iAtNum =  67 
      if(at(1:2).eq.'ER') iAtNum =  68 
      if(at(1:2).eq.'TM') iAtNum =  69 
      if(at(1:2).eq.'YB') iAtNum =  70 
      if(at(1:2).eq.'LU') iAtNum =  71 
      if(at(1:2).eq.'HF') iAtNum =  72 
      if(at(1:2).eq.'TA') iAtNum =  73 
      if(at(1:2).eq.'RE') iAtNum =  75 
      if(at(1:2).eq.'OS') iAtNum =  76 
      if(at(1:2).eq.'IR') iAtNum =  77 
      if(at(1:2).eq.'PT') iAtNum =  78 
      if(at(1:2).eq.'AU') iAtNum =  79 
      if(at(1:2).eq.'HG') iAtNum =  80 
      if(at(1:2).eq.'TL') iAtNum =  81 
      if(at(1:2).eq.'PB') iAtNum =  82 
      if(at(1:2).eq.'BI') iAtNum =  83 
      if(at(1:2).eq.'PO') iAtNum =  84 
      if(at(1:2).eq.'AT') iAtNum =  85 
      if(at(1:2).eq.'RN') iAtNum =  86 
      if(at(1:2).eq.'FR') iAtNum =  87 
      if(at(1:2).eq.'RA') iAtNum =  88 
      if(at(1:2).eq.'AC') iAtNum =  89 
      if(at(1:2).eq.'TH') iAtNum =  90 
      if(at(1:2).eq.'PA') iAtNum =  91 
      if(at(1:2).eq.'NP') iAtNum =  93 
      if(at(1:2).eq.'PU') iAtNum =  94 
      if(at(1:2).eq.'AM') iAtNum =  95 
      if(at(1:2).eq.'CM') iAtNum =  96 
      if(at(1:2).eq.'BK') iAtNum =  97 
      if(at(1:2).eq.'CF') iAtNum =  98 
      if(at(1:2).eq.'ES') iAtNum =  99 
      if(at(1:2).eq.'FM') iAtNum = 100 
      if(at(1:2).eq.'MD') iAtNum = 101 
      if(at(1:2).eq.'NO') iAtNum = 102 
      if(at(1:2).eq.'LR') iAtNum = 103 
      if(at(1:2).eq.'RF') iAtNum = 104 
      if(at(1:2).eq.'DB') iAtNum = 105 
      if(at(1:2).eq.'SG') iAtNum = 106 
      if(at(1:2).eq.'BH') iAtNum = 107 
      if(at(1:2).eq.'HS') iAtNum = 108 
      if(at(1:2).eq.'MT') iAtNum = 109 
c      if(iAtNum.eq.0) then
c        write(*,857) atom
c857     format(' Atom ',a,' is not supported. Program termiiAtNumed...')
c      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c last non-blank character in line
      Function LENSTR(string)
      implicit double precision (a-h,o-z)
      Character*(*) string
      Do 1000 i=Len(string),1,-1
         If (string(i:i) .Ne. ' ') Then
            LENSTR = i
            Return
         Endif
 1000 Continue
      LENSTR = 0
      Return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determines the last non-blank character in line
      integer function ils(line)
      character*(*) line
      ils=len(line)
      do i=1,len(line)      
        if(line(i:i).ne.' ') ils=i
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c determines the first non-blank character in line
      integer function ifs(line)
      character*(*) line
      ifs=0
      do i=1,len(line)      
        if(line(i:i).ne.' ') then
          ifs=i
          return
        endif  
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! e.g. Cu(12) --> Cu12 etc.

      subroutine erbras(atstr)

        implicit none

!------ change length probably
        character(len=10) :: atstr
        integer :: i
        
        intent(inout) :: atstr

        !write(*,*) atstr 

!------ cut "("
        do i=1,10 ! <--
          if(atstr(i:i)=='(') then
            atstr=atstr(1:i-1)//atstr(i+1:10) ! <--
            exit
          end if
        end do

!------ cut ")"
        do i=1,10 ! <--
          if(atstr(i:i)==')') then
            atstr=atstr(1:i-1)//atstr(i+1:10) ! <--
            exit
          end if
        end do

!------ just in case
        atstr=adjustl(trim(atstr))

        !write(*,*) atstr

        return

      end subroutine erbras

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
! e.g. Cu(12) --> Cu etc.

      subroutine ernum(atstr)

        implicit none

!------ change length probably
        character(len=10) :: atstr
        integer :: i
        
        intent(inout) :: atstr

        !write(*,*) atstr 

!------ cut "("
        do i=1,10 ! <--
          if(atstr(i:i)=='(') then
            atstr=atstr(1:i-1) ! <--
            exit
          end if
        end do

!------ just in case
        atstr=adjustl(trim(atstr))

        !write(*,*) atstr

        return

      end subroutine ernum

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
! e.g. Cu(12) --> Cu12 etc.

      subroutine erUM(atstr)

        implicit none

!------ change length probably
        character(len=10) :: atstr
        integer :: i
        
        intent(inout) :: atstr

        !write(*,*) atstr 

!------ cut "("
        do i=1,10 ! <--
          if(atstr(i:i+1)=='UM') then
            atstr=atstr(1:i-1)//atstr(i+2:10) ! <--
            exit
          end if
        end do

!------ just in case
        atstr=adjustl(trim(atstr))

        !write(*,*) atstr

        return

      end subroutine erUM

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      




