// v1: By Luiz Guimaraes Ferreira     At some point in the past      guima00@gmail.com    
// v2: By Filipe Matusalem            March 2014                     filipematus@gmail.com
//   Change log:
//     - Changed coding language from Fortran to C
// v3: By Bruno Lucatto               March 2017                     brunolucatto@gmail.com
//   Change log:
//     - Corrected the starting point of the summation on the fourier transform
//     - Added comments explaining most passages
//     - Minor code changes to improve readability and remove unecessary steps
//     - Small correction on the value of the constants
//
// notes:
//   - Uses the output of atm_cGuima3 program
//   - For this code, V_spinup = V_spindown
//   - INPUTS to this code must be named s follow:
//     - INP.ae-05
//     - VTOTAL.ae
//     - VTOTAL.ae-05
//     - POTCAR
//   - some tests suggests that it does not matter wether to use VTOTAL0 or VOTAL1 files,
//     as long as you keep the same number for both the .ae and .ae-05 files.
//   - Compile with:   $:g++ -o add2POTCARv3 add2POTCAR-v3.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Function to jump lines to auxiliate reading files
void jumplines(FILE *file, int nlines);

int main(int argc, char *argv[]){

  // Declaring variables
  FILE *potcar,*potcarnew,*vfile,*vfile05,*inp;
  int i,j,a,jumps,Nr,Nk,orbital,core,valence,n;
  unsigned long int pos1,pos2;
  double r[1500],Vatomo[1500],Vion[1500],Vr[1500],Vk[1500];
  double kmax,k,Ry,a0,pi,fourier_j,amplitude,occupation,CUT;
  char ch, str[200], str1[50], Linha[200];

  // Constants
  Ry=13.605803; // Rydberg constant in eV
  a0=0.52917706; // Bohr radius in Angstroms
  pi=3.141592654;
  
  // Trimming function exponent (we always use 8)
  n = 8;

  // Warning in case of incorrect call
  if( argc < 2 ){
    printf("\n\n");
    printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR \n");
    printf("EE        RR    RR   RR    RR   OO   OO   RR    RR \n");
    printf("EE        RR    RR   RR    RR   OO   OO   RR    RR \n");
    printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR \n");
    printf("EE        RRRR       RRRR       OO   OO   RRRR     \n");
    printf("EE        RR  RR     RR  RR     OO   OO   RR  RR   \n");
    printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR \n\n");

    printf("Enter the CUT value as argument for the program! e.g. ./add2POTCARv3 3.70\n"); 
    printf("The amplitude can be entered as a second argument, default A = 1.0 e.g. ./add2POTCARv3 3.70 0.8\n");
    printf("\n\n");
    exit(0);
  }

  // CUT parameter for the trimming function
  CUT = atof(argv[1]);    

  // Amplitude of the potential. Default: A = 1.0
  if(argc > 2)
    amplitude = atof(argv[2]);
  else {
    amplitude = 1;
  }
  
  //////////////////////////////////// OPENING AND CREATING FILES ////////////////////////////////////
  
  // Opens POTCAR file, read only
  potcar = fopen("POTCAR","r"); 
  if(!potcar){
    printf("Error opening POTCAR file\n");
    exit(0);
  }

  // Creates a new POTCAR file and names it properly
  if(amplitude == 1.0)
    sprintf(str, "POTCARcut%.3f",CUT);
  else  
    sprintf(str, "POTCARcut%.3fA%.1f",CUT,amplitude);
  potcarnew = fopen(str,"w"); 
  if(!potcarnew){
    printf("Error creating POTCARcut file%s\n",str);
    exit(0);
  }

  // Opens the potential of the atom
  vfile = fopen("VTOTAL.ae","r"); 
  if(!vfile){
    printf( "Error opening VTOTAL1.ae file\n");
    exit(0);
  }

  // Opens the potential of the ion
  vfile05 = fopen("VTOTAL.ae-05","r");
  if(!vfile05){
    printf( "Error opening VTOTAL1.ae-05 file\n");
    exit(0);
  }

  // Opens the input for the ionic run of the atm_cGuima3 program
  inp = fopen("INP.ae-05","r");
  if(!inp){
    printf( "Error opening INP.ae-05 file\n");
    exit(0);
  }

  //////////////////////////////////// READING THE DATA FROM THE INPUTS ////////////////////////////////////

  // atm_cGuima3 potentials --------------------------------------------------------------------------

  // Places the cursor at the beginning of the data of the spin-down potential and
  // counts how many values for the radius there are. Ps: Only spin down is ready,
  // since it is assumed that V_spinup = V_spindown.
  jumplines(vfile,1);
  Nr=0;
  fscanf(vfile,"%s",str);
  while(strcmp(str,"Down")!=0){
    fscanf(vfile,"%s",str);
    Nr++;
  }                               
  jumplines(vfile,2);
  
  // Reads the radius values of the ion's potential and places the cursor 
  // at the beginning of the data of the spin-down potential
  jumplines(vfile05,1);
  for(i=0;i<Nr;i++){ 
    fscanf(vfile05,"%lg",&r[i]);
  }
  jumplines(vfile05,3);

  // Reads both the atom's and ion's potentials
  for(i=0;i<Nr;i++){
    fscanf(vfile,"%lg",&Vatomo[i]);
    fscanf(vfile05,"%lg",&Vion[i]);
  }
  
  // Applies the trimming function to the difference between the potentials and multiplies by
  // the pre-factor, so it's Fourier transform can be directly added to the data on POTCAR.
  for(i=0;i<Nr;i++){
    if(r[i]<CUT)
      Vr[i] = 4*pi*Ry*pow(a0,3) * pow(1-pow(r[i]/CUT,n),3) * (Vion[i]-Vatomo[i]) * amplitude;
    else
      Vr[i] = 0; 
  }

  // Original POTCAR file -----------------------------------------------------------------------------
  
  // Places the cursor at the end of the description of the parameters used to generate the POTCAR,
  // stores its position and places it on the beginning of the file.
  while(strcmp(str,"PSCTR-controll")!=0)
    fscanf(potcar,"%s",str);
  pos1 = ftell(potcar);   
  rewind(potcar);
  // Counts the number of line jumps that are necessary to reach that position on the POTCAR file and
  // add two jumps to this value. Hence, by jumping this amount of lines, the cursor will be at the 
  // beginning of the "local part" data. Finally, places the cursor at the beginning of the file.
  jumps=0;
  pos2=0;
  while(pos2<=pos1){
    ch = getc(potcar);
    if(ch=='\n')
      jumps++;
    pos2 = ftell(potcar);
  }
  jumps=jumps+2;
  rewind(potcar);
  
  // Input for the ionic run of the atm_cGuima3 program -------------------------------------------------

  jumplines(inp,5);
        
  // Discards the number of core orbitals, reads the number of valence orbitals and jumps one line,
  // in case there are any comments
  fscanf(inp,"%*d %d",&valence);
  jumplines(inp,1);

  strcpy(str,""); // Erases data on str
  for(i=0;i<valence;i++){
    fscanf(inp,"%d",&core);        // Reads the energy level number of the orbital (n)
    fscanf(inp,"%d",&orbital);     // Reads the angular momentum number of the orbital (l)
    fscanf(inp,"%lg",&occupation); // Reads its occupation
    jumplines(inp,1);              // Places the cursor in the next line

    // Creates a string with the information about the occupation used to generate the correction
    if(orbital==0)sprintf(str1,"%ds%.3lg ",core,occupation);
    if(orbital==1)sprintf(str1,"%dp%.3lg ",core,occupation);
    if(orbital==2)sprintf(str1,"%dd%.3lg ",core,occupation);
    if(orbital==3)sprintf(str1,"%df%.3lg ",core,occupation);
	  strcat(str,str1);
    
  }
  /////////////////////////////////////////// WRITING THE NEW POTCAR //////////////////////////////////////////
  
  // Reads the first line of the original POTCAR and print it to the new file, adding the 
  // information about the correction to it.
  fscanf(potcar,"%200[^\n]",Linha);
  fprintf(potcarnew,"%s CUT=%.3f %sAmp. = %.1f LDA-1/2 %s",Linha,CUT,str,amplitude,__DATE__);  

  // Copies all the information until the "local part" data to the new POTCAR
  for(i=0;i<jumps;i++){
    do{
      ch = getc(potcar);
      fprintf(potcarnew,"%c",ch);}           
    while(ch!='\n');
  }
  
  // Stores the kmax value and copies it to the new POTCAR
  fscanf(potcar,"%lg",&kmax);
  fprintf(potcarnew," %.18lg\n",kmax);

  
  // Counts the number of Fourier coefficients while stores the value of the potential for each of them in a vector.
  // The first fscanf performed has Nk=0, as expected. The update is done before the fscanf so when the reading is
  // not succesfull, Nk is not incremented.
  Nk=-1;
  do{
    Nk++;
    a=fscanf(potcar,"%lg",&Vk[Nk]);}
  while(a!=0);
  
  // Calculates and prints the corrected values to the new POTCAR -----------------------------------------------------
  for(j=0;j<Nk;j++){

    if(j==0)
      k=1e-12; // Takes approximately the limit of sin(x)/x when x->0, while getting rid of the pole at k=0
    else
      k=j*kmax/Nk;
    
    // The pseudopotential is given in terms of the radial distance, and is only defined for r >= 0, as expected. Since it is
    // only evaluated inside an integral from 0 to infinity, it does not matter what values it assumes for r < 0. A natural
    // choice is to define the function to be zero for negative values, but a more convenient choice is to choose
    // v(-r)=-v(r) and n(-r)=-n(r), since purely real and odd functions have purely imaginary Fourier transforms.
    // Let v' and n' be the odd extensions of the potential and the number density, respectively. 
    //                              
    //      /\ inf          /\ inf                  /\ inf                   /\ inf
    //      |               |                       |                        |
    // Ev = |  v(r)n(r)dr = |  v'(r)n'(r)dr = (1/2) |  v'(r)n'(r)dr = -(1/2) | V(k)N(k)dk
    //      |               |                       |                        |
    //     \/  0           \/  0                   \/ -inf                  \/ -inf
    // 
    // On the third equalitty, we used the fact that the product of two odd functions is even, and in the last step
    // we have applied Parseval's theorem, considering that the Fourier transforms are purely imaginary. Even though the
    // function may not pass through the origin, we can still make an odd extension, by making it discontinuous.
    //
    // The data stored on POTCAR corresponds to the Fourier transform of the odd extension of v. It can be approximated
    // by the summation on the right, where the prefactors were ommited.
    //                                                  ______
    //                     /\ inf                      \     | Nr
    //                     |                            \
    // V(k) = i*sqrt(2/pi) |  v(r)sin(bkr)dr  => V(k) ~  >     (v[i]sin(bkr[i]) + v[i-1]sin(bkr[i-1]))/2 * (r[i]-r[i-1])
    //                     |                            /
    //                    \/ 0                         /_____| 1
    //          
    // Computes the opposite of the imaginary part of the j-th fourier transform coefficient through numerical integration
    // Index zero stands for the r=DeltaR, and the function is assumed to be zero at the origin. Thus, the first trapezium
    // of the numerical integration is degenerated to a triangulum, and its area must be calculated as so.
    
    fourier_j = Vr[0]*sin(a0*k*r[0])/2 * r[0];
    for(i=1;r[i]<CUT && i<Nr;i++){ // The second condition is added just for safety
      fourier_j += (Vr[i]*sin(a0*k*r[i])+Vr[i-1]*sin(a0*k*r[i-1]))/2 * (r[i]-r[i-1]);
    }

    Vk[j] = Vk[j] + fourier_j/(a0*k);       // Adds the two coefficients
    if(j>0 && j%5==0)fprintf(potcarnew,"\n"); // Line break each five columms
    fprintf(potcarnew," %15.8E",Vk[j]);      // Print corrected value to the new POTCAR
  }

  fprintf(potcarnew,"\n ");
  
  // Copies all other data to the new POTCAR
  while((ch = getc(potcar)) != EOF){
    fprintf(potcarnew,"%c",ch);
  }
  
  fclose(potcar);
  fclose(potcarnew);
  fclose(vfile);
  fclose(vfile05);
  fclose(inp);   
}

/////////////////////////////////////////// JUMP LINES //////////////////////////////////////////
// This function jumps a specifed amount of lines, and places the cursor at the beginning of
// the next line (i.e. after the '\n')

void jumplines(FILE *file, int nlines){
  int i,j;
  char ch;
  const int MAXlength=1000;
  
  for(i=0;i<nlines;i++){
    
    fscanf(file,"%c",&ch);
    for(j=0;ch!='\n' && ch!=EOF && j<MAXlength;j++)
      fscanf(file,"%c",&ch);
    
    if(ch==EOF){
      printf("\nERROR: 'jumplines' routine reached End Of File\n");
      exit(0);
    }
    else if(j==MAXlength){
      printf("\nERROR: Too many characters on a single line (>%d). Check the inputs.\n",MAXlength);
      exit(0);
    }
    
  } 
}