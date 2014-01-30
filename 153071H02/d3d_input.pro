PRO d3d_input,inputs                                   
user=+GETENV('USER')
;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=153071											;; Shot Number
time=3.290 		        							;; Time (s) 
runid='153071H02'     		; runid of FIDASIM
device='D3D'										;; D3D,NSTX,AUGD,MAST
result_dir='/u/stagnerl/FIDASIM/RESULTS/D3D/'  	;; Location where results will be stored result_dir/runid directory will be created
profile_dir='/u/heidbrin/FERRON/'              	;; Location of profile save files. EX: profile_dir+'shot/'+'dne142353.00505'

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/u/heidbrin/FERRON/153071/153071H01_fi_x.cdf'  ;; CDF file from transp with the distribution function
emin=5.                                                      ;; minimum energy used from the distribution function
emax=100.      ;; maximum energy used from the distribution function
pmin=-1.       ;; Minimum pitch used from the distribution function
pmax=1.        ;; Maximum pitch used from the distribution function

;;-----------------------------------------------------
;; Beam/diagnostic/equilibrium Selection
;;-----------------------------------------------------
isource=5     		    ;; Beam source index (FIDASIM only simulates one NBI source)
einj=81.3                 ;; [keV] If 0, get data from MDS+
pinj=2.56                 ;; [MW] If 0, get data from MDS+

diag=['OBLIQUE','NPA']	;; Name of the FIDA diag

gfile=''                ;; If empty, use MDS+; otherwise, filename
equil='EFIT08'			;; MDS+ equilibrium runid

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=1.d0		;; Use -1 when Bt and Ip are in the opposite direction   
ab=2.01410178d0     ;; Atomic mass of beam [u]
ai=2.01410178d0     ;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=6 	;; 5: BORON, 6: carbon, 7: Nitrogen

;;-----------------------------------------------------
;; Wavelength Grid
;;-----------------------------------------------------
lambdamin=6470.d0       						;; Minimum wavelength of wavelength grid[A] 
lambdamax=6670.d0       						;; Maximum wavelength of wavelength grid[A] 
nlambda=2000L           						;; Number of wavelengths
dlambda= (lambdamax-lambdamin)/double(nlambda)	;; Wavelength separation

;;---------------------------------------------------
;; Define FIDASIM grid in machine coordinates(x,y,z)
;;---------------------------------------------------
origin=[-141.09,-230.4,0.0]
alpha=1.3D ;; rotation about z axis
beta=0.0d
nx=60				;; Number of cells in x direction
ny=40				;; Number of cells in y direction
nz=50				;; Number of cells in z direction
xdim1=45.			;; Minimum x value
xdim2=150.			;; Maximum x value
ydim1=-40.			;; Minimum y value
ydim2=60.			;; Maximum y value
zdim1=-70.			;; Minimum z value
zdim2=70.			;; Maximum z value

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fast=5000000   			;; FIDA
nr_ndmc=50000 			;; Beam emission
nr_halo=500000	   		;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
ne_wght=50   				;; Number of Pitches, energyies and gyro angles 
np_wght=50   				;; Number of Pitches, energyies and gyro angles 
nphi_wght=50   				;; Number of Pitches, energyies and gyro angles 
emax_wght=100  				;; Maximum energy (keV)
ichan_wght=-1  				;; -1 for all channels, otherwise a given channel index
dwav_wght=.1   				;; Wavelength interval
wavel_start_wght=651.  		;; Minimum wavelength
wavel_end_wght=663.   		;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[0]   					;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[1]   ;; (0 or 1) If 1 then spectra are calculated
calc_birth=[0]  ;; (0 or 1) If 1 then the birth profile is calculated
f90brems=[0]    ;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_fida_wght=[1]   ;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[1]    ;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]    ;; (0 or 1) If 1 then the neutral density is loaded from an existing 
load_fbm=[1]         ;; (0 or 1) If 1 then the fbm is loaded (over written by calc_spec or calc_npa) 
ps=[0]          ;; (0 or 1) If 1 then make hard copy of plots

;;------------------------------------------------
;; DO NOT MODIFY THIS PART
;;------------------------------------------------
install_dir=+GETENV('FIDASIM_DIR')
inputs={shot:shot,time:time,runid:runid,device:strupcase(device),install_dir:install_dir,result_dir:result_dir,$
        cdf_file:cdf_file,profile_dir:profile_dir,emin:emin,emax:emax,pmin:pmin,pmax:pmax,isource:isource,diag:diag,$
        einj:einj,pinj:pinj,equil:equil,btipsign:btipsign,ab:ab,ai:ai,impurity_charge:impurity_charge,$
        lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,$
        nx:nx,ny:ny,nz:nz,xdim1:xdim1,xdim2:xdim2,ydim1:ydim1,ydim2:ydim2,zdim1:zdim1,zdim2:zdim2,$
        origin:origin,alpha:alpha,beta:beta,nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo,$
        ne_wght:ne_wght,np_wght:np_wght,nphi_wght:nphi_wght,$
        emax_wght:emax_wght,ichan_wght:ichan_wght,dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
        wavel_end_wght:wavel_end_wght,calc_npa:calc_npa,calc_spec:calc_spec,calc_birth:calc_birth,calc_fida_wght:calc_fida_wght,$
        calc_npa_wght:calc_npa_wght,f90brems:f90brems,load_neutrals:load_neutrals,load_fbm:load_fbm,ps:ps}


END

