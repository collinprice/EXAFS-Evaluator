#include "file_read_write.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_xyzrc(char *fname)
{
  FILE *fp;
  float xt,yt,zt,ct;
  int i;

  ct=0;
  fp = fopen(fname,"wb");
  for(i=0;i<N;i++) 
    {
      xt=X[i];yt=Y[i];zt=Z[i];
      fwrite(&xt, sizeof(float),1,fp);
      fwrite(&yt, sizeof(float),1,fp);
      fwrite(&zt, sizeof(float), 1,fp);
      fwrite(&radius[i], sizeof(float), 1,fp);
      fwrite(&ct, sizeof(float), 1,fp);
    }
  fclose(fp);
}


void center(void)
{
  long int i;
  double Xt,Yt,Zt;
  Xt=Yt=Zt=0.0;

 for(i=0;i<N;i++)
   {
     Xt+=X[i]; 
     Yt+=Y[i]; 
     Zt+=Z[i];
   }
 Xt/=N; Yt/=N; Zt/=N;
 printf("System center is at: %f %f %f\n", Xt,Yt,Zt);

 for(i=0;i<N;i++)
   {
     X[i]-=Xt; 
     Y[i]-=Yt; 
     Z[i]-=Zt; 
   }

}



void remove_bad_wat(long start, float zcen, float z_width, float margin)
{
  long i,j,c;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  
  xmin=ymin=zmin=1000;
  xmax=ymax=zmax=-1000;

 // Find minimum and maximum along each of the axes
  for(i=0;i<N;i++)
    {
      if(X[i]<xmin)xmin=X[i];
      if(Y[i]<ymin)ymin=Y[i];
      if(Z[i]<zmin)zmin=Z[i];
      if(X[i]>xmax)xmax=X[i];
      if(Y[i]>ymax)ymax=Y[i];
      if(Z[i]>zmax)zmax=Z[i];
    }
  printf("\nXmin=%.3f Xmax=%.3f\n",xmin, xmax);
  printf("Ymin=%.3f Ymax=%.3f\n",ymin, ymax);
  printf("Zmin=%.3f Zmax=%.3f\n",zmin, zmax);
  printf("Cell size:            %.2f %.2f %.2f\n",xmax-xmin,ymax-ymin,zmax-zmin); 
  xmin+=margin;ymin+=margin;xmax-=margin;ymax-=margin;
  printf("Cropped cell will be: %.2f %.2f %.2f\n",xmax-xmin,ymax-ymin,zmax-zmin); 
  
  for(j=0;j<NRES;j++)
    good_list[j]=1;

  for(j=start,c=0;j<NRES;j++)
    {
      for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
	if( ((Z[i] > zcen - z_width ) && (Z[i] < zcen + z_width)) 
	    | (X[i] < xmin ) |  (X[i] > xmax) |  (Y[i] < ymin ) | (Y[i] > ymax))
	  {
	    good_list[j]=0;
	    c++;
	    break;
	  }       
    }  
  printf("Removed %li residues\n",c);
}


void read_gaulog(char *gaulog)
{
  FILE *fp,*fp1;
  int i,j,n;
  char line_buf[90];
  char label[50];

 fp=fopen(gaulog, "rt");
  if(fp==NULL)
    {printf("** E ** gaussian log not found\n");return;}

  for(n=1;;)
    {
      if(!fgets(line_buf,82,fp)) // Read line from log file
	break; 
      
      for(i=0;i<strlen(line_buf);i++)
	if(line_buf[i]!=' ')break;
 
     if(!strncmp(&line_buf[i],"Standard orientation",20))
	{
	  sprintf(label,"set.%i.rst7",n);
	  fp1=fopen(label,"wt");
	  fgets(line_buf,82,fp);
	  fgets(line_buf,82,fp);
	  fgets(line_buf,82,fp);
	  fgets(line_buf,82,fp);  
	  
	  i=0;
	  for(;;)
	    {
	      fgets(line_buf,82,fp);
	      if(sscanf(&line_buf[35],"%lf%lf%lf",&X[i],&Y[i],&Z[i])!=3)break;
	      i++;
	    }	  
	  N=i;n++;
	 
	  fprintf(fp1,"\n %li\n ",N);
	  for(i=0,j=0;i<N;i++)
	    {
	    fprintf(fp1,"%12.7f%12.7f%12.7f", X[i],Y[i],Z[i]); j++;
	    if(j==2)
	      {fprintf(fp1,"\n ");j=0;}
	    }
	  fclose(fp1);
	}
    }
  printf("NATOMS=%li\n",N);
  printf("Extracted %i coordinate sets\n",n-1);
  fclose(fp);
}


void read_resp_charges(char *respcrg)
{
 char line_buf[200];  
 int i;
 FILE *fp;

 fp=fopen(respcrg, "rt");
  if(fp==NULL)
    {printf("** E ** topology file not found\n");return;}
  /* ----------------- READ CHARGES ---------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[10],"Point Charges",13))
	break;
    }
  fgets(line_buf,82,fp);
  fgets(line_buf,82,fp);
  for(i=0;i<N;i++)
    {
      fgets(line_buf,82,fp);
      sscanf(&line_buf[30],"%f",&charge[i]); 
    }
  fclose(fp);
}


int read_parm7(char *prmtop)
{
  FILE *fp;
  long i,j,k;
  long junk;
  char line_buf[200];  

  fp=fopen(prmtop, "rt");
  if(fp==NULL)
    {printf("** E ** topology file not found\n");return(0);}
  fgets(line_buf,82,fp);
  if(strncmp(&line_buf[1],"VERSION",7))
    {printf("** E ** topology file is not in parm7 format\n");return(0);}

  /* ----------------- READ NATOMS ---------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG POINTERS",13))
	break;
    }
  fgets(line_buf,82,fp);
  fscanf(fp,"%li%li%li%li",&N,&junk,&NBONH,&NBONA);  //read NATOM
  fgets(line_buf,82,fp);
  fscanf(fp,"%li%li",&i,&NRES);
  printf("Reading AMBER topology file:\n");
  printf("NATOM=%li, NRES=%li, NBONH=%li, NBONA=%li\n",N,NRES,NBONH,NBONA);
  /*------------------ READ ATOM NAMES -----------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG ATOM_NAME",14))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0,k=0;i<=N/20;i++)
    {
      fgets(line_buf,82,fp);
      for(j=0;j<20;j++)
	if(k<N)
	  {strncpy(atom_name[k],&line_buf[j*4],4);k++;}
	else 
	  break;
      if(k==N)break;
    }
  /*------------------- READ CHARGES -----------------*/
  while(1)
      {
	fgets(line_buf,82,fp);
	if(!strncmp(&line_buf[1],"FLAG CHARGE",11))
	  break;
      } 
  fgets(line_buf,82,fp);
  for(i=0;i<N;i++)
    fscanf(fp,"%f",&charge[i]);
  // Convert to electron charges
  for(i=0;i<N;i++)
    charge[i]/=18.222300; 
  /*------------------- READ MASS -----------------*/
  while(1)
      {
	fgets(line_buf,82,fp);
	if(!strncmp(&line_buf[1],"FLAG MASS",9))
	  break;
      }
  fgets(line_buf,82,fp);
  for(i=0;i<N;i++)
    fscanf(fp,"%f",&mass[i]); 
  /*----------------- READ RESIDUE LABELS --------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG RESIDUE_LABEL",14))
	break;
    }
  fgets(line_buf,82,fp);
  
  int halffull=0;
  if(NRES%20)halffull=1;
  for(i=0,k=0;i < (trunc(NRES/20)+halffull);i++)
    {
      fgets(line_buf,82,fp);
      for(j=0;j<20;j++)
	if(k<NRES)
	  strncpy(RESIDUE_LABEL[k++],&line_buf[j*4],3);
	else 
	  break;
    } 
 
  /*------------------ READ RESIDUE POINTERS ------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG RESIDUE_POINTER",20))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0;i<NRES;i++)
    {
    fscanf(fp,"%li",&RESIDUE_POINTER[i]);
    }
  RESIDUE_POINTER[NRES]=N+1;
  
  /*---- Convert RESUDUE_POINTERS and RESIDUE_LABELS to PDB ----*/
  for(i=0;i<NRES;i++)
    for(j=RESIDUE_POINTER[i]-1;j<RESIDUE_POINTER[i+1]-1;j++)
      {
	resSeq[j]=i;
	strncpy(resName[j],RESIDUE_LABEL[i],3);
      }
  /*---------------- READ BONDS INCLUDING HYDROGENS -------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG BONDS_INC_HYDROGEN",23))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0;i<NBONH;i++)
    {
      fscanf(fp,"%li%li%li",&IBH[i],&JBH[i],&junk);
      IBH[i]/=3;
      JBH[i]/=3;
    }
  /*---------------- READ BONDS WITHOUT HYDROGENS -------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG BONDS_WITHOUT_HYDROGEN",27))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0;i<NBONA;i++)
    {
      fscanf(fp,"%li%li%li",&IB[i],&JB[i],&junk);
      IB[i]/=3;
      JB[i]/=3; 
    }
  /*------------------ READ AMBER ATOM TYPES -------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG AMBER_ATOM_TYPE",20))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0,k=0;i<ceil(N/20);i++)
    {
      fgets(line_buf,82,fp);
      for(j=0;j<20;j++)
	if(k<N)
	  strncpy(Type[k++],&line_buf[j*4],2);
	else 
	  break;
    }
  /*------------------ READ RADII ---------------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(&line_buf[1],"FLAG RADII",10))
	break;
    }
  fgets(line_buf,82,fp);
  for(i=0;i<N;i++)
    fscanf(fp,"%f",&radius[i]);
  fclose(fp);
  return(1);
}


void write_xyzr(char *label)
{
  FILE *fp;
  long int i;

  fp=fopen(label, "w");
  if(fp==NULL)
    {printf("** E ** cannot write xyzr\n");exit(0);}
  
  for(i=0;i<N;i++)
    fprintf(fp, "%f %f %f %f\n", X[i], Y[i], Z[i],radius[i]);
  
 fclose(fp);
}

void write_mcce_ff(char *label)
{
  int i;
  FILE *fp;

  //VDW_RAD  NTR01  CA  1.908
  //VDW_EPS  NTR01  CA  0.1094



  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}  
  printf("Writing MCCE force field file: %s\n", label);

  for(i=0;i<N;i++)
    {

    fprintf(fp,"%s  %s %s %s %i\n", "VDW_RAD",  RESIDUE_LABEL[0],atom_name[i], Type[i], connect_num[i]);
    fprintf(fp,"%s  %s %s %s %i\n", "VDW_EPS",  RESIDUE_LABEL[0],atom_name[i], Type[i], connect_num[i]);
    }
  fclose(fp);

}





void write_mcce(char *label)
{
  long int i,j;
  int D_chrg;
  FILE *fp;
  float chrg, chrgf;

  chrg=0;
  chrgf=0;

  // MAKE LISTS OF BONDED ATOMS
  for(i=0;i<N;i++)
    {// MCCE reads only 3 digits in charges  
      chrg+=nearbyint(charge[i]*1000)/1000;
      chrgf+=charge[i];
      connect_num[i]=0;
      // BONDS WITHOUT HYDROGEN
      for(j=0;j<NBONA;j++)
	{
	  if(IB[j]==i)
	    connect_list[i][connect_num[i]++]=JB[j];  
	  if(JB[j]==i)
	    connect_list[i][connect_num[i]++]=IB[j];
	}
      // BONDS WITH HYDROGEN     
      for(j=0;j<NBONH;j++)
	{
	  if(IBH[j]==i)
	    connect_list[i][connect_num[i]++]=JBH[j];
	  if(JBH[j]==i)
	    connect_list[i][connect_num[i]++]=IBH[j];
	}  
    }

  D_chrg = nearbyint((chrg-chrgf)*1000);

  if(nearbyint(chrgf*1000)/1000 != 0.0)
    printf("\n ** WARNING! ** Total charge is %.4f\n Check if the charge was intended to be non zero\n\n", chrgf);
  // CORRECT CHARGE
  
  printf("Total/Trunc charge: %.4f/%.3f, Diff: 0.00%i\n",chrgf, chrg, D_chrg );
  printf("Changing charges on %i atoms by 0.001\n", D_chrg );
  for(i=0;i<abs(D_chrg);i++)
    if(D_chrg>0)
      charge[N-1-i]-=0.001;
    else
      charge[N-1-i]+=0.001;
  chrg=0;
  for(i=0;i<N;i++)
    chrg+=nearbyint(charge[i]*1000)/1000;
  printf("Corrected charge:   %.4f\n", chrg);

  wrap_names();
  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}  
  printf("Writing MCCE topology file: %s\n", label);
  
  fprintf(fp,"CONFLIST %s        %sBK\n\n",RESIDUE_LABEL[0],RESIDUE_LABEL[0]);
  fprintf(fp,"NATOM    %sBK     %3li\n\n",RESIDUE_LABEL[0],N);
  
  for(i=0;i<N;i++)
    if(strcspn(atom_name[i]," ")==4)
      fprintf(fp,"IATOM    %sBK %s %li\n",RESIDUE_LABEL[0],atom_name[i],i);
    else
      fprintf(fp,"IATOM    %sBK  %s%li\n",RESIDUE_LABEL[0],atom_name[i],i);
     
  fprintf(fp,"\n");
  
  for(i=0;i<N;i++)
    if(strcspn(atom_name[i]," ")==4)
      fprintf(fp,"ATOMNAME %sBK %4li %s\n",RESIDUE_LABEL[0],i,atom_name[i]);
    else
      fprintf(fp,"ATOMNAME %sBK %4li  %s\n",RESIDUE_LABEL[0],i,atom_name[i]);
  
  fprintf(fp, "\n#1.Basic Conformer Information: name, pka, em, rxn\n");
  fprintf(fp, "PROTON   %sBK      0 \n\n",RESIDUE_LABEL[0]);
  fprintf(fp, "PKA      %sBK      0.0 \n\n",RESIDUE_LABEL[0]);
  fprintf(fp, "ELECTRON %sBK      0 \n\n",RESIDUE_LABEL[0]);
  fprintf(fp, "EM       %sBK      0.0 \n\n",RESIDUE_LABEL[0]);
  fprintf(fp, "RXN      %sBK      0.0 \n\n",RESIDUE_LABEL[0]);
  
  fprintf(fp, "\n#2.Structure Connectivity:\n");
  // PRINT CONNECTIVITY
  for(i=0;i<N;i++)
    {
      if(strcspn(atom_name[i]," ")==4)
	fprintf(fp,"CONNECT  %sBK %s ",RESIDUE_LABEL[0],atom_name[i]);
      else 
	fprintf(fp,"CONNECT  %sBK  %s",RESIDUE_LABEL[0],atom_name[i]);
  
      if(connect_num[i]==4)
	fprintf(fp,"sp3       ");
      if((connect_num[i]==2)|(connect_num[i]==3))
	fprintf(fp,"sp2       ");
      if(connect_num[i]==1)
	fprintf(fp,"s         ");
     
      for(j=0;j<connect_num[i];j++)
	if(strcspn(atom_name[connect_list[i][j]]," ")==4)
	  fprintf(fp,"0    %s ",atom_name[connect_list[i][j]]);		
	else 
	  fprintf(fp,"0     %s",atom_name[connect_list[i][j]]);	
      fprintf(fp,"\n");
    }
 
  fprintf(fp, "\n#3.Atom Parameters: Partial Charges and Radii\n");
  for(i=0;i<N;i++)
    if(strcspn(atom_name[i]," ")==4)
      fprintf(fp,"RADIUS   %s   %s %.2f\n",RESIDUE_LABEL[0],atom_name[i],radius[i]);	     		
    else 
      fprintf(fp,"RADIUS   %s    %s%.2f\n",RESIDUE_LABEL[0],atom_name[i],radius[i]);	   
  
  fprintf(fp,"\n");
  for(i=0;i<N;i++)
    if(strcspn(atom_name[i]," ")==4)
      fprintf(fp,"CHARGE   %sBK %s %8.3f\n",RESIDUE_LABEL[0],atom_name[i],nearbyint(charge[i]*1000)/1000);	  	     		
    else 
      fprintf(fp,"CHARGE   %sBK  %s%8.3f\n",RESIDUE_LABEL[0],atom_name[i],nearbyint(charge[i]*1000)/1000);	 	   
  fclose(fp);




}





void read_msms_vert(char *label, double *x, double *y, double *z)
{
  /*------- MSMS format description ------*/
  /*------------- 3 Header lines ---------*/
  /* 1st and 2nd line:  comments          */
  /* 3rd line contain 4 numbers:          */    
  /* 1 the number of triangles,           */
  /* 2 the number of spheres in the set,  */
  /* 3 the triangulation density          */
  /* 4 the probe sphere radius.           */ 
  /* ---------------------- VERTEX DATA -----------------------*/
  /* 1-3 the coordinates (x,y,z) and the normals (nx,ny,nz)    */
  /* 4-6 the number of the face to which the vertex belongs:   */
  /* 7   0 - the vertex belongs to the analytical surface      */
  /*     negative - the vertex is on the edges of this surface */ 
  /* 8   Index of the closest sphere.                          */
  /* 9   1 for vertices which belong to toric reentrant faces  */
  /*     (including vertices of the analytical surface)        */
  /*     2 for vertices inside reentrant faces                 */
  /*     3 for vertices inside contact faces.                  */

  long i, ntri, nsph; 
  double rho, prad, nv_x, nv_y, nv_z;
  char line_buf[100];
  int surf, face_type;
  long sph_index;

  FILE *fp;
  printf("Reading msms vertices:\n");
  
  fp=fopen(label, "rt");
  if(fp==NULL)
    {printf("** E ** cannot read msms vertex\n");exit(0);}
  
  fgets(line_buf,82,fp);
  printf("%s",line_buf);
  fgets(line_buf,82,fp);
  
  fscanf(fp,"%li%li%lf%lf",&ntri, &nsph, &rho, &prad);
  printf("ntri=%li nsph=%li rho=%f prad=%f\n",ntri, nsph, rho, prad);
  
  for(i=0;i<ntri;i++)
    fscanf(fp,"%lf%lf%lf%lf%lf%lf%i%li%i",&x[i],&y[i],&z[i],&nv_x,&nv_y,&nv_z,&surf,&sph_index,&face_type);
  fclose(fp);
  
  // Write vertices in .gro format 
  printf("\nWriting vertices in .gro format\n");
  fp=fopen("vertices.gro", "wt");
  fprintf(fp,"%s\n%li\n","MSMS vertex",ntri);
  for(i=0;i<ntri;i++)
    fprintf(fp,"%5lx%5s%5s%5lx%8.3f%8.3f%8.3f\n",i,"HULL","C",i,X[i]/10,Y[i]/10,Z[i]/10);
  for(i=0;i<9;i++)
    fprintf(fp,"%10.5f",1.0);
  fprintf(fp,"\n"); 
  fclose(fp);
}

void read_pdb(char *label)
{
  unsigned long  n, i;
  char *p;
  char resID[5];
  char line_buf[200];  
  int flag_dec_limit, flag_end_of_9999_residue;

  FILE *fp;
  fp = fopen(label,"rt");
  if(fp == NULL){printf("** E ** file %s not found, abort\n", label); exit(1);}
  printf("Reading PDB file\n");
  n=0;
  flag_dec_limit=0;
  flag_end_of_9999_residue=0;
  while(1){
    if(!fgets(line_buf,80,fp)) // Read line from PDB file
      break; 
    if(!strncmp(&line_buf[0],"END",3))
      break;
    // Ignore everything except ATOM and HETATM
    if(!strncmp(&line_buf[0],"ATOM",4) | !strncmp(&line_buf[0],"HETATM",6))
      { 
	strncpy(atom[n],&line_buf[0],6);
	//		printf("%s",atom[n]);           
	sscanf(&line_buf[6],"%5li",&serial[n]);
	//		printf("%5i",serial[n]);
	strncpy(atom_name[n], &line_buf[12], 4);
	//		printf("%s",atom_name[n]);
	strncpy(altLoc[n],&line_buf[16], 1);
	//		printf("%s",altLoc[n]);
	strncpy(resName[n],&line_buf[17], 3);
	//		printf("%s",resName[n]);
	//		printf(" ");
	strncpy(chainID[n],&line_buf[21], 1);
	//		printf("%s",chainID[n]);
	strncpy(resID,&line_buf[22], 4);
	resSeq[n]=strtol(resID,&p,10);
	if(resSeq[n]==9999)
	  // DEC limit reached
	  flag_dec_limit=1; 
	if((resSeq[n]<resSeq[n-1]) && flag_dec_limit)
	  // End atom of residue
	  flag_end_of_9999_residue=1;
	// if DEC limit reached and END of residue
	if(flag_end_of_9999_residue && flag_dec_limit)
	  resSeq[n]=strtol(resID,&p,16);
	strncpy(iCode[n],&line_buf[26], 1);
	sscanf(&line_buf[30],"%lf%lf%lf%f%f",&X[n], &Y[n], &Z[n], &occupancy[n], &tempFactor[n]);
	n++;
	N=n;
      }
  }
  fclose(fp);

  resSeqN[0]=0;
  for(i=1;i<N;i++) // Make actual list of residues valid for nres>10000
    {
      resSeqN[i]=resSeqN[i-1];
      if(resSeq[i]!=resSeq[i-1])
	resSeqN[i]++;
    }
  printf("NATOM=%li atoms, NRES=%li\n", N, resSeqN[N-1]+1);
}


void read_namd_binary(char *label, unsigned long n, double *x, double *y, double *z)
{
  int i,na;
  FILE *fp;

  printf("Reading NAMD binary:\n");
  fp=fopen(label, "r");
  if(fp==NULL)
    {printf("** E ** cannot read namd binary file\n");exit(0);}
  
  fread(&na,4,1,fp); 
  printf("NATOM=%i\n",na);
  if((int) (na-n) != 0){
    printf("Number of atoms does not match topology\n");
    exit(0);}

  for(i=0;i<n;i++)
    {
      fread(&x[i],8,1,fp);
      fread(&y[i],8,1,fp);
      fread(&z[i],8,1,fp);
    }

  fclose(fp);
}


void insert_wat_namdbin(char *label, double *x, double *y, double *z)
{
  // Inserts water atoms (O,H1,H2) with coordinates x,y,z


}


void check_bonds(void)
{
  float R, B_max, BH_max;
  long j;
  int count;
 
  B_max=1.9;
  BH_max=1.55;
  printf("Bonds longer than %f/%f will be printed\n\n",B_max,BH_max);

      // BONDS WITHOUT HYDROGEN
  count=1;
      for(j=0;j<NBONA;j++)
	{
	  R=sqrt( (X[IB[j]]-X[JB[j]])*(X[IB[j]]-X[JB[j]])+
		  (Y[IB[j]]-Y[JB[j]])*(Y[IB[j]]-Y[JB[j]])+
		  (Z[IB[j]]-Z[JB[j]])*(Z[IB[j]]-Z[JB[j]]));
	  if(R>B_max)
	    printf("%i Bond %f between atoms [%li %li] %s-%li: %s-%s\n",
		   count++,R,IB[j],JB[j],resName[IB[j]],resSeq[IB[j]], atom_name[IB[j]],atom_name[JB[j]]);
    
	}
      // BONDS WITH HYDROGEN     
      for(j=0;j<NBONH;j++)
	{
	  R=sqrt( (X[IBH[j]]-X[JBH[j]])*(X[IBH[j]]-X[JBH[j]])+
		  (Y[IBH[j]]-Y[JBH[j]])*(Y[IBH[j]]-Y[JBH[j]])+
		  (Z[IBH[j]]-Z[JBH[j]])*(Z[IBH[j]]-Z[JBH[j]]));
	  if(R>BH_max)
  printf("%i Bond %f between atoms [%li %li] %s-%li: %s-%s\n",
	 count++,R,IBH[j],JBH[j],resName[IBH[j]],resSeq[IBH[j]], atom_name[IBH[j]],atom_name[JBH[j]]);	  	 
	}  
}


void add_wat(double x, double y, double z)
{
  NRES++;
  N++;
  RESIDUE_POINTER[NRES]=N+1;
  strcpy(atom_name[N-1], "O   ");
  strcpy(RESIDUE_LABEL[NRES-1],"WAT");
  strcpy(resName[N-1],RESIDUE_LABEL[NRES-1]);
  X[N-1]=x;
  Y[N-1]=y;
  Z[N-1]=z;
  good_list[NRES-1]=1;
  printf("Added WATER resid #%li xyz: %.3f %.3f %.3f   \n",NRES, x,y,z);
}


void read_amber_coor(char *label, long n, double *x, double *y, double *z)
{
  long i,na;
  char line_buf[200];
  FILE *fp;

  fp=fopen(label, "rt");
  if(fp==NULL)
    {printf("** E ** cannot read coor file\n");exit(0);}
  
  fgets(line_buf,82,fp);
  fscanf(fp,"%li",&na);

 if(n!=na){
    printf("Number of atoms does not match prmtop\n");
    exit(0);}

  for(i=0;i<n;i++)
    fscanf(fp,"%lf%lf%lf",&x[i],&y[i],&z[i]);
  
  fclose(fp);
}

void wrap_names(void)
{
  int i;
  char tmp[4];

  for(i=0;i<N;i++)
    {
      if(strcspn(atom_name[i]," ")==4)
	{
	  strncpy(&tmp[1],atom_name[i],3);
	  strncpy(atom_name[i],&atom_name[i][3],1);
	  strncpy(&atom_name[i][1],&tmp[1],3);
	}
    }
}


void convert_amb_types_to_radii(void)
{
  int i;
  for(i=0;i<N;i++)
    {      
      /* CARBON with one hydrogen or aromatic */
      if(!strncmp("ab",Type[i],2)|
	 !strncmp("sb",Type[i],2)|
	 !strncmp("bb",Type[i],2)|
	 !strncmp("cs",Type[i],2)|
	 !strncmp("tb",Type[i],2)|
	 !strncmp("pb",Type[i],2)|
	 !strncmp("qb",Type[i],2)|
	 !strncmp("rb",Type[i],2)|
	 !strncmp("zk",Type[i],2)|
	 !strncmp("za",Type[i],2)|
	 !strncmp("ze",Type[i],2)|
	 !strncmp("q2",Type[i],2)|
	 !strncmp("qq",Type[i],2)|
	 !strncmp("k1",Type[i],2)|
	 !strncmp("k2",Type[i],2)|
	 !strncmp("k3",Type[i],2)|
	 !strncmp("k4",Type[i],2)|
	 !strncmp("k5",Type[i],2)|
	 !strncmp("qo",Type[i],2)|
	 !strncmp("EC",Type[i],2)|
	 !strncmp("C ",Type[i],2)|
	 !strncmp("CA",Type[i],2)|
	 !strncmp("CB",Type[i],2)|
	 !strncmp("CC",Type[i],2)|
	 !strncmp("CD",Type[i],2)|
	 !strncmp("CK",Type[i],2)|
	 !strncmp("CM",Type[i],2)|
	 !strncmp("CN",Type[i],2)|
	 !strncmp("CQ",Type[i],2)|
	 !strncmp("CR",Type[i],2)|
	 !strncmp("CV",Type[i],2)|
	 !strncmp("CW",Type[i],2)|
	 !strncmp("C*",Type[i],2)|
	 !strncmp("CY",Type[i],2)|
	 !strncmp("CZ",Type[i],2)|
	 !strncmp("qo",Type[i],2)) 
	{radius[i]=1.7;continue;}  
      /* CARBON with 2 or three hydrogens */
      if(!strncmp("c2",Type[i],2)|
	 !strncmp("t1",Type[i],2)|
	 !strncmp("t2",Type[i],2)|
	 !strncmp("t3",Type[i],2)|
	 !strncmp("CT",Type[i],2))
	{radius[i]=2.0;continue;}
      /* SULPHUR */
      if(!strncmp("S ",Type[i],2)|
	 !strncmp("SH",Type[i],2))
	{radius[i]=2.7;continue;}
     /* PHOSPHORUS */
      if(!strncmp("P ",Type[i],2))
	{radius[i]=2.8;continue;}
      /* ALIPHATIC HYDROGEN */
      if(!strncmp("H ",Type[i],2)|
	 !strncmp("HC",Type[i],2)|
	 !strncmp("H1",Type[i],2)|
	 !strncmp("H2",Type[i],2)|
	 !strncmp("H3",Type[i],2)|
	 !strncmp("HA",Type[i],2)|
	 !strncmp("H4",Type[i],2)|
	 !strncmp("H5",Type[i],2)|
	 !strncmp("HO",Type[i],2)|
	 !strncmp("HS",Type[i],2)|
	 !strncmp("HW",Type[i],2)|
	 !strncmp("HP",Type[i],2)|
	 !strncmp("HZ",Type[i],2)|
	 !strncmp("hc",Type[i],2)|
	 !strncmp("hn",Type[i],2)|
	 !strncmp("ha",Type[i],2))
 	{radius[i]=1.2;continue;}
      /* NITROGEN */
      if(
	 !strncmp("N ",Type[i],2)|
	 !strncmp("NA",Type[i],2)|
	 !strncmp("NB",Type[i],2)|
	 !strncmp("NC",Type[i],2)|
	 !strncmp("N2",Type[i],2)|
	 !strncmp("N3",Type[i],2)|
	 !strncmp("NT",Type[i],2)|
	 !strncmp("N*",Type[i],2)|
	 !strncmp("NY",Type[i],2)|
	 !strncmp("ns",Type[i],2)|
	 !strncmp("nh",Type[i],2)|
	 !strncmp("nm",Type[i],2))
 	{radius[i]=1.55;continue;}
      /* OXYGEN */
      if(!strncmp("O ",Type[i],2)|
	 !strncmp("O2",Type[i],2)|
	 !strncmp("OG",Type[i],2)|
	 !strncmp("OW",Type[i],2)|
	 !strncmp("OH",Type[i],2)|
	 !strncmp("OS",Type[i],2)|
	 !strncmp("o1",Type[i],2)|
	 !strncmp("o2",Type[i],2)|
	 !strncmp("o3",Type[i],2)|
	 !strncmp("o4",Type[i],2))
 	{radius[i]=1.52;continue;}
      /* MAGNESIUM */
      if(!strncmp("mx",Type[i],2))
	{radius[i]=1.18;continue;}
      /* IRON */
      if(!strncmp("fn",Type[i],2)|
	 !strncmp("fh",Type[i],2) )
	{radius[i]=1.48;continue;}
   /* DUMMY */
      if(!strncmp("du",Type[i],2))
	{radius[i]=0.0;continue;}
      /* MANGANESE -? */
      if(!strncmp("mo",Type[i],2)|
	 !strncmp("m1",Type[i],2)|
	 !strncmp("m2",Type[i],2)|
	 !strncmp("m3",Type[i],2)|
	 !strncmp("m4",Type[i],2))
	{radius[i]=1.9;continue;}
      /* CALCIUM ION -? */
      if(!strncmp("ca",Type[i],2))
	{radius[i]=2.23;continue;}
      /* SODIUM ION -? */
      if(!strncmp("IP",Type[i],2))
	{radius[i]=2.0;continue;}

	printf("ERROR: Type %s not found\n",Type[i]);
	exit(0);
   
    }
}
  


void write_namd_binary(char *label, long n, double *x, double *y, double *z)
{
  long int i;
  FILE *fp;

  fp=fopen(label, "w");
  if(fp==NULL)
    {printf("** E ** cannot write namd binary file\n");exit(0);}

  fwrite(&n,4,1,fp);
  for(i=0;i<n;i++)
    {
      fwrite(&x[i],8,1,fp);
      fwrite(&y[i],8,1,fp);
      fwrite(&z[i],8,1,fp);
    }
  fclose(fp);
}



void write_one_pdb(char *label)
{
  unsigned long int i,j;
  unsigned long resCount;
  FILE *fp;
  
 
  // printf("Writing file: %s\n",label);
  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}
  resCount=0;
 
  for(j=0;j<NRES;j++)
    {	  
      for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)	  
	if(resCount>9998)
	  {
	    if(i>99998)
	      { // Print atom numbers in HEX, resid in HEX
		if(strlen(atom_name[i])<2)
		  fprintf(fp, "ATOM  %6lx%s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6lx %s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	    else
	      { // Print atom numbers in DEC, resid in HEX
		if(strlen(atom_name[i])<2)
		  fprintf(fp, "ATOM  %6li%s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6li %s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	  }
	else
	  {
	    if(i>99998)
	      { // Print atom numbers in HEX, resid in DEC
		if(strlen(atom_name[i])<2)
		  fprintf(fp, "ATOM  %6lx%s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6lx %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	    else
	      { // Print atom numbers in DEC, resid in DEC
		if(strlen(atom_name[i])<2)
		  fprintf(fp, "ATOM  %6li%s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	  }
      resCount++;
    }    
  fprintf(fp,"END\n");
  if(resCount>65534)
    printf("DUDE, the number of residues is huge!\n File will not comply with PDB standard\n"); 
  fclose(fp);
}




void write_one_pdb_sel(char *label)
{
  unsigned long int i,j;
  unsigned long resCount;
  FILE *fp;
  
 
  printf("Writing file: %s\n",label);
  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}
  resCount=0;
 
  for(j=0;j<NRES;j++)
    {
      if(good_list[j])
	{	  
	  for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)

	    if(resCount>9998)
	      {
		if(i>99998)
		  { // Print atom numbers in HEX, resid in HEX
		    if(strlen(atom_name[i])<2)
		      fprintf(fp, "ATOM  %6lx%s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		    else
		      fprintf(fp, "ATOM %6lx %s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		  }
		else
		  { // Print atom numbers in DEC, resid in HEX
		    if(strlen(atom_name[i])<2)
		      fprintf(fp, "ATOM  %6li%s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		    else
		    fprintf(fp, "ATOM %6li %s %s %5lx    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			    i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		  }
	      }
	    else
	      {
		if(i>99998)
		  { // Print atom numbers in HEX, resid in DEC
		    if(strlen(atom_name[i])<2)
		      fprintf(fp, "ATOM  %6lx%s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		    else
		      fprintf(fp, "ATOM %6lx %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		  }
		else
		  { // Print atom numbers in DEC, resid in DEC
		    if(strlen(atom_name[i])<2)
		      fprintf(fp, "ATOM  %6li%s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			      i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		    else
		    fprintf(fp, "ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			    i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		  }
	      }
	  resCount++;
	}    
    }
  fprintf(fp,"END\n");
  if(resCount>65534)
    printf("DUDE, the number of residues is huge!\n File will not comply with PDB standard\n"); 
  fclose(fp);
}



void write_split_pdb_sel(char *label)
{
  unsigned long int i,j;
  unsigned long chunk_count,resCount;
  char fullabel[80];
  FILE *fp;
  
  chunk_count=1; i=0; resCount=0;
  
  strcpy(fullabel,label);
  sprintf(&fullabel[strlen(label)],".%li",chunk_count);
  printf("Writing file: %s\n",fullabel);
  fp = fopen(fullabel,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}
  resCount=0;

  for(j=0;j<NRES;j++)
    {
      if(resCount>9999)
	{ // Start new pdb file if residue NUM > 9999
	  fprintf(fp,"END\n");
	  fclose(fp);
	  resCount=0;
	  chunk_count++;
	  sprintf(&fullabel[strlen(label)],".%li",chunk_count);
	  printf("Writing file: %s\n",fullabel);
	  fp = fopen(fullabel,"wt");
	  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}
	}
      if(good_list[j])
	{	  
	  for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
	    {
	    if(i>99999)
	      { // Print atom numbers in HEX
		if(strlen(atom_name[i])<2)
		  fprintf(fp, "ATOM  %6lx%s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6lx %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	    else
	      { // Print atom numbers in DEC
		if(strcspn(atom_name[i]," ")==4)
		  fprintf(fp, "ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
		else
		  fprintf(fp, "ATOM %6li  %s%s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
			  i+1, atom_name[i], resName[i], resCount+1, X[i], Y[i], Z[i], occupancy[i], tempFactor[i]);
	      }
	  // Add TER records to the ends of protein chains
	  if(!strncmp(atom_name[i],"OXT",3))
	    fprintf(fp,"TER\n");
	    }
	  resCount++;
	}    
    }
  fprintf(fp,"END\n");
  fclose(fp);
}

void mass_to_element(void) {
   long int i;
  for(i=0;i<N;i++)
    {
      // printf("MASS #%f\n",mass[i]);
      if((int)nearbyint(mass[i])==1)
        {strncpy(Element[i],"H",1);continue;}
      if((int)nearbyint(mass[i])==12)
        {strncpy(Element[i],"C",1);continue;}
      if((int)nearbyint(mass[i])==14)
        {strncpy(Element[i],"N",1);continue;}
      if((int)nearbyint(mass[i])==16)
        {strncpy(Element[i],"O",1);continue;}
       if((int)nearbyint(mass[i])==23)
        {strncpy(Element[i],"Na",2);continue;}
      if((int)nearbyint(mass[i])==24)
        {strncpy(Element[i],"Mg",2);continue;}
      if((int)nearbyint(mass[i])==32)
        {strncpy(Element[i],"S",1);continue;}
      if((mass[i]>54.99)&&(mass[i]<55.1))
        {strncpy(Element[i],"Fe",2);continue;}
      if((mass[i]>54.5)&&(mass[i]<55.0))
        {strncpy(Element[i],"Mn",2);continue;}
      if((int)nearbyint(mass[i])==35)
        {strncpy(Element[i],"Cl",2);continue;}
      if((int)nearbyint(mass[i])==65)
        {strncpy(Element[i],"Zn",2);continue;}
      if((int)nearbyint(mass[i])==40)
        {strncpy(Element[i],"Ca",2);continue;}

      printf("UNKNOWN ELEMENT #%li\n",i);
    }
}
