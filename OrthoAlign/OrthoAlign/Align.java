
import java.util.ArrayList;
import java.util.Arrays;


public class Align

	{

	
		    int []X; 
		    int []Y;
		    
		    String []st1;  // s�quence X originale sans les,
		    
		    String []st2; // s�quence Y originale sans les,
		     
		    
		    
		    
		    String ancetre="";
		    
		    
		    
		    
		    String al1, al2; // les deux lignes de l'alignement
		    
		    int tablePD [][]; // table de la programmation dynamique
		    
		    int infini=99999;
		    
		    int cout_cyc=0, cout_acyc=0, cout; // cout de l'alignement cyclique et acyclique
		    
		    long time=0; // temps de calcul de l'alignement
		    
		    String hhmmss; // temps de calcul sous la forme hh:mm:ss
		    
		    boolean hascycle=false;
		    
		    
		    
		    int tmp_min=infini;
		    
		    int minValue[][]; // minimum entre tous les evenements consideres
		    
		    
		    int nbMatch=0, nbLX=0, nbLY=0, nbDupX=0, nbDupXY=0, nbDupIX=0, nbDupIXY=0;
		    int nbDupY=0, nbDupYX=0, nbDupIY=0, nbDupIYX=0, nbInvXY=0, nbSub=0, nbTransp=0;
		    
		    
		    int stat_dup []=new int [100];
		    
		    int stat_loss []=new int [100];
		    
		    int stat_inver []=new int [100];
		    
		    int stat_trans []=new int [100];
		    
		    boolean dup_perteX [], dup_perteY [];
		    
		    
		    
		    Occurrence maxdupX[], maxdupY[];
		    // tableaux des plus longues duplications de X(resp. Y) avec source dans X (resp. dans Y)
		    
		    
		    Occurrence maxdupXY[], maxdupYX[]; 
		    // tableaux des plus longues duplications de X(resp. Y) avec source dans Y (resp. dans X)
		    
		    
		    Occurrence maxdupInvX[], maxdupInvY[];
		    // tableaux des maximum duplications inverses de X (resp. Y) avec source dans X (resp. dans Y) 
		    
		    
		    Occurrence maxdupInvXY[], maxdupInvYX[];
		    // tableaux des maximum duplications inverses de X (resp. Y) avec source dans Y (resp. dans X)
		    
		    
		    
		    //debut de la partie ajoutee le 24nov
		    
		    ArrayList<Integer> maxInvXY[][]; // tableaux des maximum inversions entre X et Y
		    
		  //fin de la partie ajoutee le 24nov
		    
		    
		    Occurrence minInvXY[][]; // tableaux des minimum inversions entre X et Y
		    
		    //boolean singletonX[], singletonY[];
		    
		    int ifDupX[][];  //valeur de l'alignement si (i,j) est une duplication de X avec source dans X
		    int ifDupY[][];  //valeur de l'alignement si (i,j) est une duplication de Y avec source dans Y
		    int ifLX[][];    //valeur de l'alignement si (i,j) est une perte dans X
		    int ifLY[][];    //valeur de l'alignement si (i,j) est une perte dans Y
		    int ifMatch[][]; //valeur de l'alignement si (i,j) est un match
		    
		    int ifDupXY[][]; //valeur de l'alignement si (i,j) est une duplication de X avec source dans Y
		    int ifDupYX[][]; //valeur de l'alignement si (i,j) est une duplication de Y avec source dans X
		    
		    int ifDupInvX[][]; // valeur de l'alignement si (i,j) est une duplication
		                       // inverse de X avec source dans X
		    
		    int ifDupInvY[][]; // valeur de l'alignement si (i,j) est une duplication
		    				   // inverse de Y avec source dans Y
		    
		    
		    int ifDupInvXY[][]; // valeur de l'alignement si (i,j) est une duplication
			   					// inverse de X avec source dans Y
		    
		    
		    int ifDupInvYX[][]; // valeur de l'alignement si (i,j) est une duplication
		    					// inverse de Y avec source dans X
		    
		    
		  //debut de la partie ajoutee le 24nov
		    
		    int ifInvXY[][]; // valeur de l'alignement si (i,j) est une inversion
		    
		    
		  //fin de la partie ajoutee le 24nov
		    
		    
		    int posDupX[][];     // position de la duplication dans X qui donne le meilleur alignement 
		    					// exemple : si le cout de ifDX[i][j]=C[x][j]+c(D(i-x)) avec X[l..i]\
		    					// est la plus plongue duplication qui commence par i alors, posDupX[i][j]=x;
		    
		    
		    int posDupY[][];     // position de la duplication dans Y qui donne le meilleur alignement 
		    
		    int posLX[][];       // position de la perte dans X qui donne le meilleur alignement
		     
		    int posLY[][];       // position de la perte dans Y qui donne le meilleur alignement
		    
		    int posMatch[][];   // position du match qui donne le meilleur alignement
		    
		    int posDupXY[][];   // position de la duplication dans X (source dans Y) 
		    					//qui donne le meilleur alignement

		    
		    int posDupYX[][];   // position de la duplication dans Y (source dans X) 
								//qui donne le meilleur alignement
		    
		    
		    int posDupInvX[][]; 	// position de la duplication inverse dans X 
									//qui donne le meilleur alignement
		    
		    
		    
		    int posDupInvY[][];		// position de la duplication inverse dans Y 
									//qui donne le meilleur alignement
		    
		    
		    int posDupInvXY[][];    // position de la duplication inverse dans X (source dans Y)
		    						//qui donne le meilleur alignement
		    
		    int posDupInvYX[][];	// position de la duplication inverse dans Y (source dans X)
									//qui donne le meilleur alignement
		    
		    
		  //debut de la partie ajoutee le 24nov
		    
		    int posInvXY[][];  // taille du segment de la meilleure inversion entre X et Y
		    
		    
		  //fin de la partie ajoutee le 24nov
		    
		    
       //debut de la partie ajoutee le 28nov
		    
		    int ifSub[][]; //valeur de l'alignement si (i,j) est substitution
		    
		    
		  //fin de la partie ajoutee le 28nov
 
 
        //debut de la partie ajoutee le 28nov
		    
		    int posSub[][];	// position de la substitution qui donne le meilleur alignement
		    
		    
		     //fin de la partie ajoutee le 28nov
		      
 

		      int [] positionX;
			  
			  int [] positionY;
		    
		    
			  
			  Cover2 []isCoveredX;
			  Cover2 []isCoveredY;
			  
		    
		    
		    /*
		    
		    ArrayList<Coordinates> allDupsX; // contient toutes les duplications d'un alignement dans X
		    								// .value : taille de la duplication, .x: exptremite droite, .y: extremite gauche 	
		    
		    
		    ArrayList<Coordinates> allDupsY; // contient toutes les duplications d'un alignement dans Y
		    								 // meme commentaire	
		    
		    
		    */
		    
		    //ArrayList<Cover> covAllDupsX; // contient les couvertures de toutes les duplications d'un alignement dans X
			                                   // .value : taille de la couverture, .x: exptremite droite, .y: extremite gauche 	


		    //ArrayList<Cover> covAllDupsY; // contient les couvertures de toutes les duplications d'un alignement dans Y
			                                   // meme commentaire
		    
		    //Cover dupCoverX[], dupCoverY[];
		    
		    Cover2 dupCover[];
		    
		    
		    
		   public Align(){} ;
		    
		   public Align (int []G1, int []G2)
		   {
			   
			   this.X=G1;
			   this.Y=G2;
			   
			   
		   }
		
		   
		   
		   public Align (int []G1, int []G2, String []s1, String []s2)
		   {
			   
			   this.X=G1;
			   this.Y=G2;
			   
			   this.st1=s1;
			   this.st2=s2;
			   
			   
		   }
		
		   
		   public void set_X(int []st)
		   {
			   this.X=st;
			   
			   
		   }
		   
		   
		   public void set_Y(int []st)
		   {
			   
			   this.Y=st;
			   
		   }
		   
		   public void set_XY(int []st1, int []st2)
		   {
			   this.X=st1;
			   this.Y=st2;
			   
		   }
		   
		   
		   public int[] get_X()
		   {
			   return this.X;
			   
		   }
		   
		   
		   public int[] get_Y()
		   {
			   return this.Y;
			   
		   }
		   
		   
		   
		   
		   static int maxDup(int []X, int i) // determine la position de la plus grande duplication entre la position i et 0 sur toute la chaine de caracteres.
		    
		   {
			  //declaration des variables
		    	
		   int j, a, b, l=0, r=0, s; 
		   
		   int strLength=X.length;
		   
		 //commencer par le cote droite
		   	   
		   j=i+2; //inutile de commencer a la position i+1 
	 
		   if (i>=strLength) return -1;
		   
		   
		   while (j<strLength && r<i) //inutile d'avancer le j si le [0..i] est deja une duplication
			   
		   {
			   
			   
			   while ( j<strLength && (X[i]!=X[j]) ) 
				   
			   {
			       j++; // avancer jusqu'a egalite de caractere
			       
			   }
			       
			       
			   if ( j<strLength && (X[i-r]==X[j-r]) )
			   {
				   l=0;
				   a=i-1;
				   b=j-1;
				   
				  while(a>=0 && b>i && (X[a]==X[b])) 
				 
				    {
					  
					  l++;
					  a--;
				  	  b--;
				    }
				  
				  if (l>r) r=l;
				  
			   }
			   
			   j++;
		  }
		    	
		  
		   j=1;
		   
		   s=0;
		   
		  
		 //tester le cote gauche
		   
		   
		   
		   while (j<i-r && r<(i/2)+1) //empecher le suffixe de depasser le moitie du segment gauche+gestion des bords
			                        
		   {
			   
			   while ( j<i-r && (X[i]!=X[j]) ) 
				   
			   	{
			       j++; // reculer jusqu'a egalite de caractere
			       
			   	}
			       
			      
			   if ( (j<i-r) && (j>=r) && (X[i-r]==X[j-r]) )
			   {
				   l=0;
				   a=i-1;
				   b=j-1;
				   
				  while(a>j && b>=0 && (X[a]==X[b])) 
				 
				    {
					      
					  l++;
					  a--;
				  	  b--;
				    }
				  
				  if (l>r) r=l;
				  
			   }
			   
			   j++;
		  }
		
		   return i-r;
		      
		   }
		   
		    
		   
		   
		   
		    
		   /**
			 *
			 *  isSingleton: Verifie si le caractere X[i] est un singleton
			 */
		   
		   
		   static boolean isSingleton(int []S, int i)
			
			{
			
			   int j=(i+1)%S.length;
			   int c=S[i];
			   
			   
			   while (S[j]!=c) 
			           j=(j+1)%S.length;
			
				
				
				return (j==i);
			
			}
		  
		   
		 /**
		  * 
		  * ifMatch: retourner le score en supposant X[i]=Y[j] 
		  */
		 
		 
		   void ifMatch(int i, int j)
		 
		 	{
			 if (X[i]==Y[j] && positionX[i]==positionY[j])
				 
			 {
				 ifMatch[j+1][i+1]=tablePD[j][i];
				 
				 if (tablePD[j][i]<tmp_min)
					  
					  tmp_min=tablePD[j][i];
					 
			 }	 
				
			 else 
				 
			 	 ifMatch[j+1][i+1]=infini;
			 
		 	}
		 
		   
		   
		 //debut de la partie ajoutee le 28nov
		    
		    
		   /**
			  * 
			  * ifSubs: retourner le score si X[i]!=Y[j] 
			  */
			 
			 
		   void ifSubs(int i, int j)
		 
		 	{
			   
			   
			 if (X[i-1]!=Y[j-1] && X[i-1]!=91140 && Y[j-1]!=91140 && positionX[i-1]==positionY[j-1] )
				 
			 {
				 ifSub[j][i]=tablePD[j-1][i-1]+1;  // on consid�re le cout d'une substitution 1
				 
				 //System.out.println(ifSub[j+1][i+1]);
				 
				 if (ifSub[j][i]<tmp_min)
					  
					  tmp_min=ifSub[j][i];
					 
			 }	 
				
			 else 
				 
			 	 ifSub[j][i]=infini;
			 
		 	}
	    
	    
	     //fin de la partie ajoutee le 28nov
		   
		   
		   
		   
		 
		 /**
		  * 
		  * ifLy: retrouner le score si X[i] est une perte dans Y
		  */
		 
		   void ifLy(int i, int j)
		 
		 	{
			  
			  int min=tablePD[j][i-1]+1; // cLost(1)=1;
			  
			  ifLY[j][i]=min;
			  posLY[j][i]=i-1;
			  
			  
			  if (min<tmp_min)
			  	  tmp_min=min;
				
		 	}
		 
		 
		 /**
		  * 
		  * ifLx: retourner le score si Y[j] est une perte dans X
		  */	
		 
		 
		   
		   void ifLx(int i, int j)
		    
		    {
		    	
		    	int min=tablePD[j-1][i]+1; // cLost(1)=1;
		    	
		    	ifLX[j][i]=min;
				posLX[j][i]=j-1;
		    	
		    	
				  if (min<tmp_min)
				 
					  tmp_min=min;
				  
		    }
		  
		 
		    /**
			  * 
			  * ifDx: retourner le score si X[i] est une duplication dans X
			  */
		 
		 
		    void ifDx(int i, int j)
		  	
		  	{
			  
			  if (maxdupX[i-1].i!=-1 && X[i-1]!=91140) //if (!singletonX[i-1]) 
			  {
				  
			  int l=maxdupX[i-1].i; // car la table contient une ligne/colonne en plus
			  
			  int min=tablePD[j][l]+1;  // en realite : min=tablePD[j][l]+cDup(i-l), Dup(i-l)=1
			  
			  int x=l;
			  
			  for (int k=l+1; k<i; k++)
				  
			  	{
				  
				  int var=tablePD[j][k]+1; // +cDup(i-k);
				  
				  if (var <=min) // favoriser la plus grande chaine de duplication, on peut inverser, en mettant if (min>var) pour garder la plus grande
					             // car notre idee est de faire priorite aux matchs
					  
				  { 
					  min=var;
				      x=k;
					  
				  } 	  
					  
			  	}
				  
			  ifDupX[j][i]=min;
			  posDupX[j][i]=x;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				 
			  }
			  
			  
			  
			  }
			  
			  else ifDupX[j][i]=infini;
			  
		  	} 
		   
		    
		    
		    /**
			  * 
			  * ifDIx: retourner le score si X[i] est une duplication inverse dans X
			  */
		 
		 
		    void ifDIx(int i, int j)
		  	
		  	{
			  
			  if (maxdupInvX[i-1].i!=-1 && X[i-1]!=91140) //if (!singletonX[i-1]) 
			  {
				  
			  int l=maxdupInvX[i-1].i; // car la table contient une ligne/colonne en plus
			  
			  int min=tablePD[j][l]+1;  // en realite : min=tablePD[j][l]+cDup(i-l), Dup(i-l)=1
			  
			  int x=l;
			  
			  for (int k=l+1; k<i; k++)
				  
			  	{
				  
				  int var=tablePD[j][k]+1; // +cDup(i-k);
				  
				  if (var <=min) // favoriser la plus grande chaine de duplication, on peut inverser, en mettant if (min>var) pour garder la plus grande
					             // car notre idee est de faire priorite aux matchs
					  
				  { 
					  min=var;
				      x=k;
					  
				  } 	  
					  
			  	}
				  
			  ifDupInvX[j][i]=min;
			  posDupInvX[j][i]=x;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				 
			  }
			  
			  
			  
			  }
			  
			  else ifDupInvX[j][i]=infini;
			  
		  	} 
		    
		    
		    
		    
		    
		    /**
			  * 
			  * ifDxY: retourner le score si X[i] est une duplication dans X
			  */
		 
		 
		    void ifDxY(int i, int j)
		  	
		  	{
			  
			  if (maxdupXY[i-1].i!=-1 && X[i-1]!=91140) //if (!singletonX[i-1]) 
			  {
				  
			  int l=maxdupXY[i-1].i; // car la table contient une ligne/colonne en plus
			  
			  int min=tablePD[j][l]+1;  // en realite : min=tablePD[j][l]+cDup(i-l), Dup(i-l)=1
			  
			  int x=l;
			  
			  for (int k=l+1; k<i; k++)
				  
			  	{
				  
				  int var=tablePD[j][k]+1; // +cDup(i-k);
				  
				  if (var <=min) // favoriser la plus grande chaine de duplication, on peut inverser, en mettant if (min>var) pour garder la plus grande
					             // car notre idee est de faire priorite aux matchs
					  
				  { 
					  
					  
					  
					  min=var;
				      x=k;
					  
				  } 	  
					  
			  	}
				  
			  ifDupXY[j][i]=min;
			  posDupXY[j][i]=x;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				 
			  }
			  
			  
			  
			  }
			  
			  else ifDupXY[j][i]=infini;
			  
		  	} 
		    
		    
		    /**
			  * 
			  * ifDxY: retourner le score si X[i] est une duplication dans X
			  */
		 
		 
		    void ifInversionXY(int i, int j)
		  	
		  	{
			  
			  if (maxInvXY[j-1][i-1].size()!=0)  
			  {
			
				  
			  	  
				  
			  int l=maxInvXY[j-1][i-1].get(0); // car la table contient une ligne/colonne en plus
			  
			  int min=tablePD[j-l][i-l]+1;  // en realite : min=tablePD[j][l]+cDup(i-l), Dup(i-l)=1
			  
			  int x=l;
			  
			  for (int k:maxInvXY[j-1][i-1])
				  
			  	{
				  
				  int var=tablePD[j-k][i-k]+1; // +cInv(i-k);
				  
				  if (var <=min) 
					             
					  
				  { 
					  min=var;
				      x=k;
					  
				  } 	  
					  
			  	}
				  
			  ifInvXY[j][i]=min;
			  posInvXY[j][i]=x;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				 
			  }
			  
			  
			  
			  }
			  
			  else ifInvXY[j][i]=infini;
			  
		  	} 
		    
		    
		    
		    
		    
		    /**
			  * 
			  * ifDIxY: retourner le score si X[i] est une duplication inverse dans X avec source dans Y
			  */
		 
		 
		    void ifDIxY(int i, int j)
		  	
		  	{
			  
			  if (maxdupInvXY[i-1].i!=-1 && X[i-1]!=91140) //if (!singletonX[i-1]) 
			  {
				  
			  int l=maxdupInvXY[i-1].i; // car la table contient une ligne/colonne en plus
			  
			  int min=tablePD[j][l]+1;  // en realite : min=tablePD[j][l]+cDup(i-l), Dup(i-l)=1
			  
			  int x=l;
			  
			  for (int k=l+1; k<i; k++)
				  
			  	{
				  
				  int var=tablePD[j][k]+1; // +cDup(i-k);
				  
				  if (var <=min) // favoriser la plus grande chaine de duplication, on peut inverser, en mettant if (min>var) pour garder la plus grande
					             // car notre idee est de faire priorite aux matchs
					  
				  { 
					  min=var;
				      x=k;
					  
				  } 	  
					  
			  	}
				  
			  ifDupInvXY[j][i]=min;
			  posDupInvXY[j][i]=x;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				 
			  }
			  
			  
			  
			  }
			  
			  else ifDupInvXY[j][i]=infini;
			  
		  	}
		    
		    
		    
		  
		  /**
			  * 
			  * ifDy: retourner le score si Y[j] est une duplication dans Y
			  */
		 
		 
		  void ifDy(int i, int j)
		  	
		  	{
			  
			  if (maxdupY[j-1].i!=-1)// if (!singletonY[j-1])  
				  
			  {		  
				  
			  int l=maxdupY[j-1].i; 
			  
			  int min=tablePD[l][i]+1;// en realite : min=talblePD[l][i]+cDup(j-l), cDup(j-l)=1;
			  
			  int y=l;
			  
			  for (int k=l+1; k<j; k++)
				  
			  	{
				  
				  int var=tablePD[k][i]+1; // +cDup(j-k);
				  
				  if (var <= min) 
					  
				  {
					  min=var;
					  y=k;
					
				  }
					  
			  	}
				
			  
			  ifDupY[j][i]=min;
			  posDupY[j][i]=y;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				  
			  }
			  
			  
			  }
			  
			  else ifDupY[j][i]=infini;  
			  
		  	} 
		  
		  
		  
		  /**
			  * 
			  * ifDIy: retourner le score si Y[j] est une duplication inverse dans Y
			  */
		 
		 
		  void ifDIy(int i, int j)
		  	
		  	{
			  
			  if (maxdupInvY[j-1].i!=-1)// if (!singletonY[j-1])  
				  
			  {		  
				  
			  int l=maxdupInvY[j-1].i; 
			  
			  int min=tablePD[l][i]+1;// en realite : min=talblePD[l][i]+cDup(j-l), cDup(j-l)=1;
			  
			  int y=l;
			  
			  for (int k=l+1; k<j; k++)
				  
			  	{
				  
				  int var=tablePD[k][i]+1; // +cDup(j-k);
				  
				  if (var <= min) 
					  
				  {
					  min=var;
					  y=k;
					
				  }
					  
			  	}
				
			  
			  ifDupInvY[j][i]=min;
			  posDupInvY[j][i]=y;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				  
			  }
			  
			  
			  }
			  
			  else ifDupInvY[j][i]=infini;  
			  
		  	} 
		  
		  
		  
		  
		  /**
			  * 
			  * ifDyX: retourner le score si Y[j] est une duplication dans X
			  */
		 
		 
		  void ifDyX(int i, int j)
		  	
		  	{
			  
			  if (maxdupYX[j-1].i!=-1)// if (!singletonY[j-1])  
				  
			  {		  
				  
			  int l=maxdupYX[j-1].i; 
			  
			  int min=tablePD[l][i]+1;// en realite : min=talblePD[l][i]+cDup(j-l), cDup(j-l)=1;
			  
			  int y=l;
			  
			  for (int k=l+1; k<j; k++)
				  
			  	{
				  
				  int var=tablePD[k][i]+1; // +cDup(j-k);
				  
				  if (var <= min) 
					  
				  {
					  min=var;
					  y=k;
					
				  }
					  
			  	}
				
			  
			  ifDupYX[j][i]=min;
			  posDupYX[j][i]=y;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				  
			  }
			  
			  
			  }
			  
			  else ifDupYX[j][i]=infini;  
			  
		  	} 
		  
		  
		  /**
			  * 
			  * ifDIyX: retourner le score si Y[j] est une duplication dans X
			  */
		 
		 
		  void ifDIyX(int i, int j)
		  	
		  	{
			  
			  if (maxdupInvYX[j-1].i!=-1)// if (!singletonY[j-1])  
				  
			  {		  
				  
			  int l=maxdupInvYX[j-1].i; 
			  
			  int min=tablePD[l][i]+1;// en realite : min=talblePD[l][i]+cDup(j-l), cDup(j-l)=1;
			  
			  int y=l;
			  
			  for (int k=l+1; k<j; k++)
				  
			  	{
				  
				  int var=tablePD[k][i]+1; // +cDup(j-k);
				  
				  if (var <= min) 
					  
				  {
					  min=var;
					  y=k;
					
				  }
					  
			  	}
				
			  
			  ifDupInvYX[j][i]=min;
			  posDupInvYX[j][i]=y;
			  
			  
			  if (min<tmp_min)
			  {
				  tmp_min=min;
				  
			  }
			  
			  
			  }
			  
			  else ifDupInvYX[j][i]=infini;  
			  
		  	} 
		  
		  
		  
		  
		  
		  public static Coordinates Couverture(int []X, int d, int g)
			
			{
			    int dc=-1, gc=-1; // les deux extremites de la chaine couverture initialisees a -1
			    
			    boolean found=false, bound=false;
			    
			
				int k=X.length-1; // k prend la derniere position de la chaine 
				
				if (g>k || d<0)  // verifie s'il n y a pas de debordement 
					
					{ 
					
					    System.out.println("Message from Couverture: coordinates out of bound");
					 	return null;
					}
				
				
				while (!found && !bound)
					
				{	
				
				while( (k>g) && (k-g > g-d) && (X[g]!=X[k]) )  k--;
				
				
				if ( (k>g) && (k-g > g-d) )
					
					{
					   
					   int a=g;
					   int b=k;
					   
					   while((a>=d) && (X[a]==X[b]))
					   	{
						   
						   a--; b--;
						   
					   	}
						    
					if (a<d) // alors chaine ou couverture trouvee
						{
						       found=true; 
						       
						       dc=k-g+d;  
						       gc=k;
						    	   
						} 
					
					    
					
					else k--;
					
					}
				
				else bound=true;
				
				}
				
				
				k=d-1;
				
				while (!found && bound)
				
					{
					
					while( (k >= g-d) && (X[g]!=X[k]) )  k--;
					
					
					if ( k >= g-d )
						
					{
					
						int a=g;
						int b=k;
					   
					   while((a>=d) && (X[a]==X[b]))
					   	{
						   
						   a--; b--;
						   
					   	}
						    
					if (a<d) 
						
					{
					       found=true; 
					       
					       dc=k-g+d; 
					       gc=k;
					    	   
					} 
					  
					     else k--;
					
					}
				
					
					else bound=false; // pour sortir de la boucle
					
					
					
					}
					
				
				//System.out.println("dc= "+dc+", gc= "+gc);
			
				if (found) return new Coordinates(g-d+1, dc, gc);
				
				return null;
				
				
			
			
			}
		  
		  
		  
		  
		  /*
		   *  Chevauch(Coordinates c1, Coordinates c2) : permet de savoir s'il existe un chevauchement entre c1 et c2
		   *  
		   *  @ les coordonnees du segment du chevauchement
		   * 
		   */
		  
		  
		  
		  
		  
		  
		  public static Coordinates Chevauch(Coordinates c1, Coordinates c2)
		  
		  	{
			
			  int x1, y1, x2, y2;
			  
			  x1=c1.x;
			  y1=c1.y; // x1, y1 : coordonnees du premier segment
			  
			  
			  
			  x2=c2.x;
			  y2=c2.y; // x2, y2 : coordonnees du deuxieme segment
			  
			  
			  // une duplication de taille 1 ne peux pas etre chevauchee
			  
			  if ( (x1==y1) || (x2==y2) ) return new Coordinates(0, -1, -1); 
			  
			  
			  
			  if(y1>=y2)
				  
			  	{
				  if (y1<=x2) 
				  	{
					
					  if (x2<=x1) return new Coordinates(x2-y1+1, y1, x2);
					  else return new Coordinates(x1-y1+1, y1, x1);
					  
				  	}
					  
					  
				  else return new Coordinates(0, -1, -1);
				  
				  
				  
				  
			  	}
			  
			  else
				  
			  	{
				  if (x1>=y2) 
				  	{
					
					  if (x1<=x2) return new Coordinates(x1-y2+1, y2, x1);
					  else return new Coordinates(x2-y2+1, y2, x2);
					  
				  	}
					  
					  
				  else return new Coordinates(0, -1, -1);
				  
				  
				  
				  
			  	}
			  
		  	}
		  
		  
		  
		  public static Coordinates2 Chevauch2(Coordinates2 c1, Coordinates2 c2)
		  
		  	{
			
			  
			  
			  if ( c1.flag!=c2.flag ) return new Coordinates2(0, -1, -1,-1); 
			  
			  int x1, y1, x2, y2;
			  
			  x1=c1.x;
			  y1=c1.y; // x1, y1 : coordonnees du premier segment
			  
			  
			  
			  x2=c2.x;
			  y2=c2.y; // x2, y2 : coordonnees du deuxieme segment
			  
			  
			  // une duplication de taille 1 ne peux pas etre chevauchee
			  
			  //if ( (x1==y1) || (x2==y2) ) return new Coordinates2(0, -1, -1,-1);
			  
			  
			  if(y1>=y2)
				  
			  	{
				  if (y1<=x2) 
				  	{
					
					  if (x2<=x1) return new Coordinates2(x2-y1+1, y1, x2,c1.flag);
					  else return new Coordinates2(x1-y1+1, y1, x1,c1.flag);
					  
				  	}
					  
					  
				  else return new Coordinates2(0, -1, -1,-1); // pas de chevauchement
				  
				  
				  
				  
			  	}
			  
			  else
				  
			  	{
				  if (x1>=y2) 
				  	{
					
					  if (x1<=x2) return new Coordinates2(x1-y2+1, y2, x1,c1.flag);
					  else return new Coordinates2(x2-y2+1, y2, x2,c1.flag);
					  
				  	}
					  
					  
				  else return new Coordinates2(0, -1, -1,-1); // pas de chevauchement
				  
				  
				  
				  
			  	}
			  
		  	}
		  
		  
		  
		  /**
			 *  preKmp (x): determine le plus long prefix de [i..n] qui est suffixe de [i..n] pour chaque position dans x.
			 *
			 **/
		   
		   public static int[] preKmp4(int[] x, int p) {
			   
				int m=p+1+1;
				int kmpNext[]=new int[m];
				int i, j;

				   i = 0;
				   j = kmpNext[0] = -1;
				   
				   while (i < m-2 ) {
				      while (j > -1 && x[m-2-i] != x[m-2-j])
				         j = kmpNext[j];
				      i++;
				      j++;
				      if (x[m-2-i] == x[m-2-j])
				         kmpNext[i] = kmpNext[j];
				      else
				         kmpNext[i] = j;
				   }
				
				   
				   
				   while (j > -1 && x[m-2-i] != x[m-2-j])
				         j = kmpNext[j];
				   
				   kmpNext[m-1] = j+1;
				   
				   
					   //for (int h=m-1; h>=0; h--)
						   //System.out.print(kmpNext[h]+" ");
				 
				   
				   return kmpNext;
				   
				      
				}
		   
		   

		  
		   
		   
		   
		   
		   
		   /* createTables : Cr�ation de la table de programmation dynamique 
			   * 
			   *  Partie � developper
			   * 
			   * 
			   * 
			   * 
			   * 
			   */
		   
		   
		   
		   
		   
		  public void createTables()
		  	
		  {
			  
			  tablePD=new int[Y.length+1][X.length+1];
			  
			  int []iX=Pretraitement.inverser(X);
			  int []iY=Pretraitement.inverser(Y);
			  
			  
			  maxdupX=Pretraitement.allmaxDup(X);
			  maxdupY=Pretraitement.allmaxDup(Y);
			  
			  
			  maxdupInvX=Pretraitement.allMaxInv(X,iX);
			  maxdupInvY=Pretraitement.allMaxInv(Y,iY);
			  
			  maxdupXY=Pretraitement.allmaxDup(X,Y);
			  maxdupYX=Pretraitement.allmaxDup(Y,X);
			  
			  
			  maxdupInvXY=Pretraitement.allMaxInvP(X,iY);
			  maxdupInvYX=Pretraitement.allMaxInvP(Y,iX);
			  
			  
			 // debut de la partie ajoutee le 24nov 
			  
			  maxInvXY=Pretraitement.allLongSufInvXY2(iX, Y);
			  
			  positionX=Pretraitement.Position(X);
			  
			  positionY=Pretraitement.Position(Y);
			  
			  
			// fin de la partie ajoutee le 24nov
			  
			  
			   
			  
			  
			  
			  //singletonX=new boolean[X.length];
			  
			  //singletonY=new boolean[Y.length];
			  
			  tablePD[0][0]=0;
			  
			  
			    ifDupX=new int[Y.length+1][X.length+1];
			    ifDupY=new int[Y.length+1][X.length+1];
			    ifLX=new int[Y.length+1][X.length+1];
			    ifLY=new int[Y.length+1][X.length+1];
			    ifMatch=new int[Y.length+1][X.length+1];
			    
			    posDupX=new int[Y.length+1][X.length+1];
			    posDupY=new int[Y.length+1][X.length+1];
			    posLX=new int[Y.length+1][X.length+1];
			    posLY=new int[Y.length+1][X.length+1];
			    posMatch=new int[Y.length+1][X.length+1];
			  
			  
			  
			    //debut de la partie ajout� le 29 octobre 2013
				  
			    ifDupInvX=new int[Y.length+1][X.length+1];
			    ifDupInvY=new int[Y.length+1][X.length+1];  
			    
			    ifDupXY=new int[Y.length+1][X.length+1];
			    ifDupYX=new int[Y.length+1][X.length+1];
			    
			    ifDupInvXY=new int[Y.length+1][X.length+1];
			    ifDupInvYX=new int[Y.length+1][X.length+1];
			    
			    
			 // debut de la partie ajoutee le 24nov 
			    
			    ifInvXY=new int[Y.length+1][X.length+1];
			    
			 // fin de la partie ajoutee le 24nov 
			    
			    
				
			    posDupInvX=new int[Y.length+1][X.length+1];
			    posDupInvY=new int[Y.length+1][X.length+1];
			    
			    posDupXY=new int[Y.length+1][X.length+1];
			    posDupYX=new int[Y.length+1][X.length+1];
			    
			    
			    posDupInvXY=new int[Y.length+1][X.length+1];
			    posDupInvYX=new int[Y.length+1][X.length+1];
			    
				  
				 //fin de la partie ajout� le 29 octobre 2013
			    
			    
			    // debut de la partie ajoutee le 24nov 
			    
			    posInvXY=new int[Y.length+1][X.length+1];
			    
			   // fin de la partie ajoutee le 24nov
			    
			    
			    // debut de la partie ajoutee le 29nov 
			    
			    ifSub=new int[Y.length+1][X.length+1];
			    
			  //posSub=new int[Y.length+1][X.length+1];
				
			 // fin de la partie ajoutee le 29nov
			    
			    
			    
			    
			    
			    
			  for (int i=1; i<X.length+1; i++)
			  	{
		            
				  //maxdupXX[i-1]=maxDup(X, i-1);
				  //singletonX[i-1]=(maxdupXX[i-1]==-1);
				  
				  ifLy(i,0);
				  ifDx(i,0);
				  
				  
				  //debut de la partie ajout� le 29 octobre 2013
				  
				  ifDIx(i,0);
				  
				  ifDxY(i,0);
				  ifDIxY(i,0);
				  
				 //fin de la partie ajout� le 29 octobre 2013
				  
				  
				  
			  	tablePD[0][i]=tmp_min;
			  		
			  	
			  	tmp_min=infini;
			  	
			  	}
			 
			  
			  for (int j=1; j<Y.length+1; j++)
			  	{
				
				
				  //maxdupYY[j-1]=maxDup(Y, j-1);
				  //singletonY[j-1]=(maxdupYY[j-1]==-1);
				  
				  ifLx(0,j);
			  	  ifDy(0,j);
			  	  
			  	  
			  	  
			  	  //debut de la partie ajout� le 29 octobre 2013
				  
				  ifDIy(0,j);
				  
				  ifDyX(0,j);
				  ifDIyX(0,j);
				  
				 //fin de la partie ajout� le 29 octobre 2013
			  	  
			  	  
			  	  
			  	  
			  	  
			  	tablePD[j][0]=tmp_min;
			  		
			  	
			  	tmp_min=infini;
			  	 
			  	
			  	}
			  
			  
			  
			  for (int j=1; j<Y.length+1; j++)
				  
				  
			  {
				  
				  for (int i=1; i<X.length+1; i++)
					   
				  {  
					  
					  ifMatch(i-1,j-1);
					  
					  ifSubs(i,j);
					  
					  
					  ifInversionXY(i,j);
					  
					  ifLy(i,j);
					  ifDx(i,j);
					  ifLx(i,j);
					  ifDy(i,j);
					  
					  
					  
					//debut de la partie ajout� le 29 octobre 2013
					  
					  ifDIx(i,j);
					  
					  ifDxY(i,j);
					  ifDIxY(i,j);
					  
					  ifDIy(i,j);
					  
					  ifDyX(i,j);
					  ifDIyX(i,j);
					  
					//fin de la partie ajout� le 29 octobre 2013  
					  
					  
					  
				  	  tablePD[j][i]=tmp_min;
				  	  
				  		
				  	  tmp_min=infini;
				  	
				  }
				  
			  }
			  
			  
			  
			  
		  
		  }
		  
		  
		  
		  
		  
		  /* createAlign : Determiner un alignement possible entre X et Y, ainsi que la couverture de chaque caractere. 
		   * 
		   *  Variables : s1 : contient la premiere ligne de l'alignement, s2 : la deuxieme ligne.
		   *              s3 : contient toutes les duplications, pertes et alignement. 0: est un Match, iiii...i: est une duplication, i: perte
		   *              s4, s5 : comme s3 mais pour X (et Y) consideree seule.
		   *              k : compteur permettant de calculer la valeur de l'alignement. Ainsi, chaque duplication est de forme iiii .
		   *              kx, ky : meme principe que k mais pour X( et Y) consideree seule 
		   *              allDupsX, allDupsY : contiennent respectivement toutes les duplications dans X et Y
		   *              cpt_al : compte le nombre de caracteres de l'alignement : valeur finale : cpt+1
		   * 		      covAllDupsX, covAllDupsX ; contiennent respectivement toutes les duplications dans X(Y) et  leurs couvertures
		   *              
		   * 
		   * 
		   * 
		   * 
		   * 
		   */
		  
		  
		  public void createAlign()
		  
		  	{
			 
			  time = System.currentTimeMillis();
			  
			  createTables();
			  
			  
			  String s1="", s2="", s3="", ancetre=""; 
			  
			  /*
			   * s1    :
			   * s2   :
			   * s3   :
			   * ancetre : 
			   * 
			   */
			 
			       	
			  //String s4="", s5="" ;
			  
			  
			  int i=X.length; // repr�sente le g�nome X (colonnes de la table de programmation dynamique)
			  
			  int j=Y.length; // repr�sente le g�nome Y (lignes de la table de programmation dynamique)
			  
			  int k=0, kx=0, ky=0;  
			  
			  cout=tablePD[Y.length][X.length];
			  
			  //System.out.println("cout = "+cout);
			  
			  
			  /*
			   * k    :
			   * kx   :
			   * ky   :
			   * cout : le cout de l'alignement d�termin� par la table PD (possiblement cyclique)
			   * 
			   */
			  
			  
			  cout_cyc=cout;  // on pose cout_cyc=cout 
			  
			  /*
			  
			  ArrayList<Coordinates> allDupsX=new ArrayList<Coordinates>(); 
			  ArrayList<Coordinates> allDupsY=new ArrayList<Coordinates>();
			  
			  */
			  
			  ArrayList<Cover2> covAllDups=new ArrayList<Cover2>();
			 
			  
			  //ArrayList<Cover> covAllDupsY=new ArrayList<Cover>();
			  
			  
			  
			  //int nDup=0; //compteur du nombre de duplications dans l'alignement
			  
			  int cpt_al=-1; 
			  
			  
			  isCoveredX=new Cover2[X.length];
			  
			  isCoveredY=new Cover2[Y.length];
			  
			  dup_perteX=new boolean [isCoveredX.length];
			  dup_perteY=new boolean [isCoveredY.length];
			  
			  Arrays.fill(dup_perteX, true);
			  Arrays.fill(dup_perteY, true);
			  
			  
			  
			  
			  while (i>0 && j>0)
				  
			  {
				  
				  
				   if (ifMatch[j][i]==tablePD[j][i])  // Match
					  
				  {
					  
					  //System.out.println("match");
					  
					  cpt_al++;
					  
					 
					  int c=X[i-1];
					  
					  s1=" "+c+s1;
					  s2=" "+c+s2;
					  
					  s3="0"+s3;
					  
					  ancetre=" "+c+ancetre;
					  
					  //s4="0"+s4;
					  
					  //s5="0"+s5;
						
					  
					  isCoveredX[i-1]=new Cover2(new Coordinates2(0, i-1, i-1,0), new Coordinates2(0, j-1, j-1,1));
					  
					  isCoveredY[j-1]=new Cover2(new Coordinates2(0, j-1, j-1,1), new Coordinates2(0, i-1, i-1,0));
					  
					  i--;
					  j--;
					  
					  
					  nbMatch++;
					 
				  }
				  
				  
				// d�but de la partie ajoute�� le 24 nov
				  
				  else if (ifInvXY[j][i]==tablePD[j][i])  // Inversion
					  
				  {
					  
					  // verifie le bug si on remplace hh par posInvXY[j][i]
					  
					  //System.out.println("Je suis la");
						
					  //System.out.println("posInvXY[j][i]= "+posInvXY[j][i]);
					  
					  //System.out.println("i et j avant "+i+" ,"+j);
					  
					  int hh=posInvXY[j][i];
					  
					  int y=i;
					  int z=j;
					  
					  while (y>i-(hh))
					  
					  {	
						  
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(3, i-(hh), i-1,0), new Coordinates2(3, j-(hh), j-1,1));
						  
						  isCoveredY[z-1]=new Cover2(new Coordinates2(3, j-(hh), j-1,1), new Coordinates2(3, i-(hh), i-1,0));
						  
					 
					  y--;
					  z--;
					  
					  }
					  
					  
					  
					  j=j-(hh);
					  
					  i=i-(hh);
					  
					  
					 
					  //System.out.println("j= "+j);
					 
					  //System.out.println("i et j apres "+i+" ,"+j);
					 
					  nbInvXY++;
					  
					  
				  }
				  
				  
				// fin de la partie ajoute�� le 24 nov
				  
				  
				  
				   
				   
				  
				 // d�but de la partie ajoute�� le 23 nov
				  
				  
				  else if (ifDupXY[j][i]==tablePD[j][i])  // DXY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupXY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupXY[i-1].j-i+h+1, maxdupXY[i-1].j,1)));
					  
					  
					  //System.out.println("DUPXY--- Debut : "+h+", fin : "+(i-1));
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(7, h, i-1,0), new Coordinates2(7, maxdupXY[i-1].j-i+h+1, maxdupXY[i-1].j,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  
					  i=h;
					  
					  
					  nbDupXY++;
					  
					  
					  //System.out.println("nbDupXY : "+nbDupXY);
					  
					  
				  }
				  
				  
				  
				  
				  
				  
				  else if (ifDupInvXY[j][i]==tablePD[j][i])  // DupInvXY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvXY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupInvXY[i-1].j, maxdupInvXY[i-1].j+i-1-h,1)));
					  
					  
					  
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(8, h, i-1,0), new Coordinates2(8, maxdupInvXY[i-1].j, maxdupInvXY[i-1].j+i-1-h,1));
						  
						  
						  
						  y--;
						  
						  
					  }
					  
					  
					  i=h;
					  
					  nbDupIXY++;
					  
					  
					  
					  
				  }
				  
				  
				  
				  
				  
				 //fin de la partie ajoute�� le 23 nov
				  
				  
				  else if (ifLY[j][i]==tablePD[j][i])  // Ly
					  
				  {
					  
					  
					  
					  cpt_al++;
					  
					  //cout++;
					  
					  s1=X[i-1]+s1;
					  
					  int t=X[i-1];
					  
					  do {
						  s2="-"+s2;
						  t/=10;
						}
					  
					  
					  while (t!=0) ;
					  	
					      
					  
					  s1=" "+s1;
					  s2=" "+s2;
					  

						  
					 isCoveredX[i-1]=new Cover2(new Coordinates2(4, i-1, i-1,0), new Coordinates2(4, i-1, i-1,0));
						 
					  
					  
					  i--;
					  
					  k++;
					  
					  kx ++;
					  
					  
					  s3=Integer.toString(k)+s3;
					  
					  //s4=Integer.toString(kx)+s4;
					  
					  ancetre=" "+t+ancetre;
					  
					  nbLY++;
					  
					  
					  
				  }
				  
				  
				  
				  else if (ifDupX[j][i]==tablePD[j][i])  // Dx
					  
				  {
					  
					
					  int h=posDupX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupX[i-1].j-i+h+1, maxdupX[i-1].j,0)));
					  
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(5, h, i-1,0), new Coordinates2(5, maxdupX[i-1].j-i+h+1, maxdupX[i-1].j,0));
						  
						  y--;
						  
						  
					  }
					  
					  //cout++;
					  
					  k++;
					  
					  kx ++;
					   
					  while (i>h)
					  	
					  	{
						  cpt_al++;
						  
						  
						  s1=X[i-1]+s1;
						  
						  int t=X[i-1];
						  
						  do
						  	{
							  s2="-"+s2;
							  t/=10;
						  	}
						  
						  while (t!=0);
						  
						  
						  s1=" "+s1;
						  s2=" "+s2;
						  
						  
						  i--;
						  
						  s3=Integer.toString(k)+s3;
						  
						  //s4=Integer.toString(kx)+s4;
						  
						  
					  	}
					  
					  
					  nbDupX++;
					  
					  
					  //nDup++;
					  
				  } 
				  
				  
				  
				  else if (ifDupInvX[j][i]==tablePD[j][i])  // DupInvX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupInvX[i-1].j, maxdupInvX[i-1].j+i-h-1,0)));
					  
					  
					  
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(6, h, i-1,0), new Coordinates2(6, maxdupInvX[i-1].j, maxdupInvX[i-1].j+i-h-1,0));
						  
						  y--;
						  
						  
					  }
					  
					  
					  
					  i=h;
					  
					  
					  nbDupIX++;
					  
					  
				  }
					  
					
				    
				  
				  
				// d�but de la partie ajoute�� le 23 nov
				  
				  
				  
				  else if (ifDupYX[j][i]==tablePD[j][i])  // DYX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupYX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupYX[j-1].j-j+h+1, maxdupYX[j-1].j,0)));
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(12, h, j-1,1), new Coordinates2(12, maxdupYX[j-1].j-j+h+1, maxdupYX[j-1].j,0));
						  
						  y--;
						  
						  
					  }
					  
					  
					  j=h;
					  
					  nbDupYX++;
					  
					  
					  
				  }
				  
				  
				  
				  
				  
				  
				  
				  else if (ifDupInvYX[j][i]==tablePD[j][i])  // DupInvYX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvYX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupInvYX[j-1].j,maxdupInvYX[j-1].j+j-h-1,0)));
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(13, h, j-1,1), new Coordinates2(13, maxdupInvYX[j-1].j,maxdupInvYX[j-1].j+j-h-1,0));
						  
						  y--;
						  
						  
					  }
					  
					  
					  j=h;
					  
					  
					  nbDupIYX++;
					  
					  
				  }
				  
				  
				//fin de la partie ajoute�� le 23 nov
				  
				  
				  else if (ifLX[j][i]==tablePD[j][i])  // Lx
					  
				  {
					  
					  
					  
					  //cout++;
					  
					  cpt_al++;
					  
					  
					  int t=Y[j-1];
					  
					  do
					  	{
						  s1="-"+s1;
						  t/=10;
					  	}
					  
					  while (t!=0);
					  
					  s2=Y[j-1]+s2;
					  
					  s1=" "+s1;
					  s2=" "+s2;
					  
					  isCoveredY[j-1]=new Cover2(new Coordinates2(9, j-1, j-1,1), new Coordinates2(9, j-1, j-1,1));
					  
					  j--;
					  
					  
					  k++;
					  
					  ky ++;
					  
					  s3=Integer.toString(k)+s3;
					  
					  //s5=Integer.toString(ky)+s5;
					  
					  ancetre=" "+t+ancetre;
					  
					  
					  
					  nbLX++;
					  
					  
				  }
				  
				  
				  else  if (ifDupY[j][i]==tablePD[j][i]) // Dy
					  
				  {
					  
					  //cout++;
					  
					  int h=posDupY[j][i];
					  
					  //allDupsY.add(new Coordinates(j-h, h, j-1)); pour essayer quelque chose
					  
					  //covAllDupsY.add(new Cover(new Coordinates(j-h, h, j-1), Couverture(Y, h, j-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1),new Coordinates2(j-h, maxdupY[j-1].j-j+h+1, maxdupY[j-1].j,1)));
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(10, h, j-1,1), new Coordinates2(10, maxdupY[j-1].j-j+h+1, maxdupY[j-1].j,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  k++;
					  
					  ky ++;
					  
					  while (j>h)
					  	
					  	{
						  cpt_al++;
						  
						  int t=Y[j-1];
						  
						  do
						  	{
							  s1="-"+s1;
							  t/=10;
						  	}
						  
						  while (t!=0);
						  
						  
						  s2=Y[j-1]+s2;
						  
						  s1=" "+s1;
						  s2=" "+s2;
						  
						  
						  
						  j--;
						  
						  s3=Integer.toString(k)+s3;
						  
						  //s5=Integer.toString(ky)+s5;
						  
					  	}
					  
					  //nDup++;
					  
					  
					  nbDupY++;
					  
					  
				  }	
				  
				  
				  else if (ifDupInvY[j][i]==tablePD[j][i])  // DupInvY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupInvY[j-1].j,maxdupInvY[j-1].j+j-h-1,1)));
					  
					  
					  
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(11, h, j-1,1), new Coordinates2(11, maxdupInvY[j-1].j,maxdupInvY[j-1].j+j-h-1,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  j=h;
					  
					  nbDupIY++;
					  
				  }
				  
				  
				  // d�but de la partie ajoute�� le 29 nov
				  
				  else if (ifSub[j][i]==tablePD[j][i])  // Substitution
					  
				  {
					 	  
					  
					  //System.out.println(st1[i-1]+st2[j-1]);
					  
					  //System.out.println("Substitution entre  "+(i-1)+"et "+(j-1));
					  
					  isCoveredX[i-1]=new Cover2(new Coordinates2(1, i-1, i-1,0), new Coordinates2(1, j-1, j-1,1));
					  
					  isCoveredY[j-1]=new Cover2(new Coordinates2(1, j-1, j-1,1), new Coordinates2(1, i-1, i-1,0));
					  
					  i--;
					  j--;
					  
					  nbSub++;
					 
					  
					  
					  
				  }
				  
				  
				// fin de la partie ajoute�� le 29 nov
				  
				  
			  }
			  
			  
			  
			while (i>0)
				
				{
					
				
					
				
				//debut de la partie ajoute�� le 23 nov
				
				
				 if (ifDupXY[j][i]==tablePD[j][i])  // DXY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupXY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupXY[i-1].j-i+h+1, maxdupXY[i-1].j,1)));
					  
					  //System.out.println("DUPXY--- Debut : "+h+", fin : "+(i-1));
					  
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(7, h, i-1,0), new Coordinates2(7, maxdupXY[i-1].j-i+h+1, maxdupXY[i-1].j,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  i=h;
					  
					  nbDupXY++;
					  
					  
					  //System.out.println("nbDupXY : "+nbDupXY);
					  
					  
				  }
				  
				 
				  
				  
				  
				  else if (ifDupInvXY[j][i]==tablePD[j][i])  // DupInvXY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvXY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupInvXY[i-1].j, maxdupInvXY[i-1].j+i-1-h,1)));
					  
					  
                    int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(8, h, i-1,0), new Coordinates2(8, maxdupInvXY[i-1].j, maxdupInvXY[i-1].j+i-1-h,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  i=h;
					  
					  nbDupIXY++;
					  
					  
				  }
				
				
				
				//fin de la partie ajoute�� le 23 nov
				
				
				
				
				  else if (ifLY[j][i]==tablePD[j][i])  // Ly
					  
				  {
					
					
					//cout++;
					
					cpt_al++;
					
					
					  
					  s1=X[i-1]+s1;
					  
					  int t=X[i-1];
					  
					  do
					  	{
						  s2="-"+s2;
						  t/=10;
					  	}
					 
					  
					  while (t!=0);
					  
					  
					  s1=" "+s1;
					  s2=" "+s2;
					  
					  isCoveredX[i-1]=new Cover2(new Coordinates2(4, i-1, i-1,0), new Coordinates2(4, i-1, i-1,0));
						
					  
					  
					  i--;
					  
					  
					  k++;
					  
                  kx ++;
					  
					  
					  s3=Integer.toString(k)+s3;
					  
					  //s4=Integer.toString(kx)+s4;
					  
					  ancetre=" "+t+ancetre;
					  
					  nbLY++;
					  
				  }
				
				
				 
				 
				  else if (ifDupX[j][i]==tablePD[j][i])   // Dx
						
						
					{     
						 //cout++;
						
						  int h=posDupX[j][i];
						  
						  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour tester quelque chose
						  
						  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
						  
						  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupX[i-1].j-i+h+1, maxdupX[i-1].j,0)));
						  
						  
						  int y=i;
						  
						  while (y>h)
							  
						  {
							  
							  isCoveredX[y-1]=new Cover2(new Coordinates2(5, h, i-1,0), new Coordinates2(5, maxdupX[i-1].j-i+h+1, maxdupX[i-1].j,0));
							  
							  y--;
							  
							  
						  }
						  
						  k++;
						 
						  kx ++;
						  
						  while (i>h)
						  	
						  	{
							  cpt_al++;
							  
							  
							  
							  s1=X[i-1]+s1;
							  
							  int t=X[i-1];
							  
							  do
							  	{
								  s2="-"+s2;
								  t/=10;
							  	}
							  
							  
							  while (t!=0);
							  
							  s1=" "+s1;
							  s2=" "+s2;
							  
							  i--;
							  
							  s3=Integer.toString(k)+s3;
							  
							  //s4=Integer.toString(kx)+s4;
							  
						  	}
						  
						  
						  //nDup++;
						  
						
						  nbDupX++;
						  
					  } 
				  
				  else if (ifDupInvX[j][i]==tablePD[j][i])  // DupInvX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(i-h, h, i-1,0), new Coordinates2(i-h, maxdupInvX[i-1].j, maxdupInvX[i-1].j+i-h-1,0)));
					  
					  
					  int y=i;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredX[y-1]=new Cover2(new Coordinates2(6, h, i-1,0), new Coordinates2(6, maxdupInvX[i-1].j, maxdupInvX[i-1].j+i-h-1,0));
						  
						  
						  y--;
						  
						  
					  }
					  
					  
					  i=h;
					  
					  
					  nbDupIX++;
					  
					  
				  }
				
				
				}
			  
			  
			
			  
			while (j>0)
				
				{
				
					
					
				
				
				//debut de la partie ajoute�� le 23 nov
				
				
				 if (ifDupYX[j][i]==tablePD[j][i])  // DYX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupYX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupYX[j-1].j-j+h+1, maxdupYX[j-1].j,0)));
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(12, h, j-1,1), new Coordinates2(12, maxdupYX[j-1].j-j+h+1, maxdupYX[j-1].j,0));
						  
						  y--;
						  
						  
					  }
					  
					  
					  j=h;
					  
					  
					  nbDupYX++;
					  
					  
				  }
				  
				 
				  
				  
				  else if (ifDupInvYX[j][i]==tablePD[j][i])  // DupInvYX
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvYX[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupInvYX[j-1].j,maxdupInvYX[j-1].j+j-h-1,0)));
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(13, h, j-1,1), new Coordinates2(13, maxdupInvYX[j-1].j,maxdupInvYX[j-1].j+j-h-1,0));
						  
						  y--;
						  
						  
					  }
					  
					  j=h;
					  
					  
					  nbDupIYX++;
					  
				  }
				
				
				  else if (ifLX[j][i]==tablePD[j][i])  // Lx
					  
				  {
					
					//System.out.println("perte dans X");
					
					cpt_al++;  
					
					//cout++;
					
					int t=Y[j-1];
					  
					  do
					  	{
						  s1="-"+s1;
						  t/=10;
					  	}
					
					  while (t!=0);
					  
					  
					  s2=Y[j-1]+s2;
					  
					  s1=" "+s1;  
					  s2=" "+s2;
					  
					  
					  isCoveredY[j-1]=new Cover2(new Coordinates2(9, j-1, j-1,1), new Coordinates2(9, j-1, j-1,1));
					  
					  
					  j--;
					  
					  k++;
					  
					  ky ++;
					  
					  s3=Integer.toString(k)+s3;
					  
					  //s5=Integer.toString(ky)+s5;
					  
					  
					  ancetre=" "+t+ancetre;
					  
					  
					  
					  nbLX++;
					  
				  }
				
				
				
				
				
				//fin de la partie ajoute�� le 23 nov
				
				
				 
				  else if (ifDupY[j][i]==tablePD[j][i])    // Dy
						
						
					{
						
						//cout++;
						
						  int h=posDupY[j][i];
						  
						  //allDupsY.add(new Coordinates(j-h, h, j-1));pour essayer quelque chose
						  
						  //covAllDupsY.add(new Cover(new Coordinates(j-h, h, j-1), Couverture(Y, h, j-1)));
						  
						  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1),new Coordinates2(j-h, maxdupY[j-1].j-j+h+1, maxdupY[j-1].j,1)));
						  
						  
						  int y=j;
						  
						  while (y>h)
							  
						  {
							  
							  isCoveredY[y-1]=new Cover2(new Coordinates2(10, h, j-1,1), new Coordinates2(10, maxdupY[j-1].j-j+h+1, maxdupY[j-1].j,1));
							  
							  y--;
							  
							  
						  }
						  
						  k++;
						  
						  ky ++;
						  
						  while (j>h)
						  	
						  	{
							  
							  
							  cpt_al++;
							  
							  int t=Y[j-1];
							  
							  do
							  	{
								  s1="-"+s1;
								  t/=10;
							  	}
							  
							  while (t!=0);
							  
							  s2=Y[j-1]+s2;
							  
							  
							  s1=" "+s1;  
							  s2=" "+s2;
							  
							  
							  j--;
							  
							  s3=Integer.toString(k)+s3;
							  
							  //s5=Integer.toString(ky)+s5;
							  
							  
						  	}
						  
						  //nDup++;
						  
						 
						  nbDupY++;
						  
					  }	
				 
				 
				 
				  else if (ifDupInvY[j][i]==tablePD[j][i])  // DupInvY
					  
					  
				  {
					  //System.out.println("hahaha");
					  
					  int h=posDupInvY[j][i];
					  
					  //allDupsX.add(new Coordinates(i-h, h, i-1)); pour essayer quelque chose 
					  
					  // maxdupXX[i-1]
					  
					  //covAllDupsX.add(new Cover(new Coordinates(i-h, h, i-1), Couverture(X, h, i-1)));
					  
					  covAllDups.add(new Cover2(new Coordinates2(j-h, h, j-1,1), new Coordinates2(j-h, maxdupInvY[j-1].j,maxdupInvY[j-1].j+j-h-1,1)));
					  
					  
					  
					  
					  int y=j;
					  
					  while (y>h)
						  
					  {
						  
						  isCoveredY[y-1]=new Cover2(new Coordinates2(11, h, j-1,1), new Coordinates2(11, maxdupInvY[j-1].j,maxdupInvY[j-1].j+j-h-1,1));
						  
						  y--;
						  
						  
					  }
					  
					  
					  j=h;
					  
					  nbDupIY++;
					  
				  }
				 
					
				}
			
			
			// commentaire 1
			
			
			
			
			// Jusque la, on a un alignement et toutes les duplications possibles.
			// determinons maintenant la couverture "possiblement cyclique" de chaque duplication
			
			
			al1=s1;
			al2=s2;
			
			// masquer l'alignement pour le nouveau projet, instruction suivante
			
			// System.out.println("========Affichage de l'alignement cyclique========");
			
			//System.out.println(s3);
			
			// masquer l'alignement pour le nouveau projet, les 2 instructions suivantes
			
			//System.out.println(s1);
			//System.out.println(s2);
			
			//System.out.println(X);
			//System.out.println(s4);
			
			//System.out.println(Y);
			//System.out.println(s5);
			
			
			//System.out.println("===cout de l'alignement cyclique=== :"+cout);
			
			
			// Commentaire 2
			
			
			int nbDup=covAllDups.size(); // nbDup : le nombre de duplications
			
			
			//int nbDupY=covAllDupsY.size();
			
			
			//System.out.println("Nombre de LX avant la resolution des cycle = "+nbLX);
			
			
			
			if (nbDup!=0)
	
			{
				
				//System.out.println("Nombre de duplication dans X est : "+nbDup);
				
			dupCover=new Cover2[nbDup]; // contient toutes les duplications et leurs couvertures
			
			
			for (int var=0; var<nbDup; var++)
				
			{
					dupCover[var]=new Cover2();
					
					dupCover[var]=covAllDups.get(var);
					
					
			}
			
			
			
			// Commentaire 3
			
			
			
			
			
			//System.out.println("Affichage de toutes les duplications et sa couverture"); 
			
			
			
			//for (Cover cor:covAllDupsX)
			
			/*
			
			
			for (int var=0; var<nbDup; var++)
				
				{
						System.out.println("Duplication "+"j= "+dupCover[var].dup.y+", i= "+dupCover[var].dup.x+", Genome "+dupCover[var].dup.flag+" ,Couverture "+"j= "+dupCover[var].covDup.y+", i= "+dupCover[var].covDup.x+" ,Genome "+dupCover[var].covDup.flag);
				
				}
			
			
			*/
			
			
			
			Coordinates2 chev[][]=new Coordinates2[nbDup][nbDup]; // contient le segment en commun entre les duplications et les Cover (origines)
			
			for (int v1=0; v1<nbDup; v1++) 
				chev[v1][v1]=null;
			
			
			
			// Commentaire 4
			
			
			
			
			
			
			
			GraphAdjMat gam=new GraphAdjMat(2*nbDup, true); // le nombre de sommets est le double du nombre de duplications
			
			
			
			for (int var=0; var<nbDup; var++)
				
			{
			
				gam.insertVertex(2);
				  
				gam.insertEdge(2*var, 2*var+1); // l'origine est repr�sent� par un sommet pair, et la copie par un sommet impair
				
				  for (int va=0; va<var; va++)
				  	
				  	{
					  
					  chev[var][va]=Chevauch2(dupCover[var].dup, dupCover[va].covDup);
					  
					  chev[va][var]=Chevauch2(dupCover[va].dup, dupCover[var].covDup);
					  
					  if (chev[va][var].getValue()!=0) gam.insertEdge(2*va+1, 2*var);
					  
					  if (chev[var][va].getValue()!=0) 
						  
					  	{
						  gam.insertEdge(2*var+1, 2*va);					
						    
						    // tester le cycle
						    
						    DirectedCycle dc=new DirectedCycle(gam);
							  
							  if (dc.hasCycle())
							  	{
								  		
								  if (!hascycle) hascycle=true;
								  
								  //System.out.println("Cycle Detecte");
								  
								  ArrayList<Integer> cyc=new ArrayList<Integer>();
								  				
								   for (int c : dc.cycle())
									   
								   {
									   //System.out.print(c+" ");
									   cyc.add(c);	
								   }			
								   
								   
								   //System.out.println();
								   
								   int ctaille=cyc.size(); // taille du cycle
								   				
								   int v1=0, t=99999, sauv=-1, temp; // sauv : meilleure duplication
								   				
								   while ( (v1<ctaille-1) && (cyc.get(v1)% 2 == 0) ) v1++;
								   
								   // v1++ pour choisir un sommet impair car le chevauchement se passe entre
								   // un sommet impair et le sommet suivant dans l'ordre du cycle
								   															
								   //System.out.println("Taille du cycle = "+ctaille);
								   
								   sauv=v1; // sauvegarder le sommet impair dans sauv
								   
								   t=chev[cyc.get(v1)/2][cyc.get(v1+1)/2].value; //valeur du chevauchement possiblement =0;
								   
								   v1++; // incrementer v1 pour todo
								   
								   while (v1< ctaille-1)     // chercher un autre sommet impair et le stocker dans v2
							   		{			
								   		int v2=cyc.get(v1); // v2: sommet de duplication impair
								   						
								   		if (v2 % 2!=0)
								   			{	
								   			     
								   				temp=chev[v2/2][cyc.get(v1+1)/2].value; // valeur du chevauchement
								   				
								   				if (temp<t)
								   					
								   					{		// todo
								   					
								   					      if ((temp<t-1) || ( dupCover[v2/2].dup.value < dupCover[cyc.get(sauv)/2].dup.value) || ( dupCover[v2/2].dup.x > dupCover[cyc.get(sauv)/2].dup.x))
								   					    	  
								   					   t=temp;
								   					      
								   					   sauv=v1; // sauvegarder le meilleur sommet de duplication ayant le meilleur chevauchement avec le sommet prochain
								   					   
								   					}
								   				
								   				
								   				if (temp==t)
								   					
								   					// if ( (temp==t) && ( dupCoverX[v2/2].dup.value == temp) )
								   					
								   					if ( ( dupCover[v2/2].dup.value < dupCover[cyc.get(sauv)/2].dup.value) || ( dupCover[v2/2].dup.x > dupCover[cyc.get(sauv)/2].dup.x) )
								   					
								   					
							   					{
							   					   t=temp;
							   					   sauv=v1; 
							   					}
								   				
								   				
								   			}	
								   		
								   		v1++;
								   
							   		}
								   
								   
								   
								   int v2=cyc.get(sauv); // numero de sommet d'une duplication
								   
								   
								   // debut de la partie ajoutee le 2 dec
								   
								   Coordinates2 coorChev=chev[v2/2][cyc.get(sauv+1)/2]; // coorChev contient la partie chevauchante
								   
								   
								   int flag=coorChev.flag;
								   
								   Cover2 cov3=new Cover2();
								   
								   Cover2 cov1=new Cover2();
								   Cover2 cov2=new Cover2();
								   
								   
								   // nous avons besoin de cov1 qq la situation
								   
								   
								   Coordinates2 dup3=new Coordinates2(dupCover[v2/2].dup.value,dupCover[v2/2].dup.y,dupCover[v2/2].dup.x,dupCover[v2/2].dup.flag);	   
								    
								   Coordinates2 covdup3=new Coordinates2(dupCover[v2/2].covDup.value,dupCover[v2/2].covDup.y,dupCover[v2/2].covDup.x,dupCover[v2/2].covDup.flag);	   
									   
								   cov3.setCovDup(dup3,covdup3); // recuperer la premiere duplication
								   
								   
								   if (ctaille==5)
									   
								   {
									   
									   int a1=0;
									   
									   while ( (a1<ctaille-1) && (cyc.get(a1)% 2 == 0) ) a1++;
									   
									   
									   Coordinates2 dup1=new Coordinates2(dupCover[cyc.get(a1)/2].dup.value,dupCover[cyc.get(a1)/2].dup.y,dupCover[cyc.get(a1)/2].dup.x,dupCover[cyc.get(a1)/2].dup.flag);	   
									    
									   Coordinates2 covdup1=new Coordinates2(dupCover[cyc.get(a1)/2].covDup.value,dupCover[cyc.get(a1)/2].covDup.y,dupCover[cyc.get(a1)/2].covDup.x,dupCover[cyc.get(a1)/2].covDup.flag);	   
										   
									   cov1.setCovDup(dup1,covdup1); // recuperer la premiere duplication
									   
									   
									   // cov2 est utile seulement dans le cas des cycles de taille 5
									     
									   Coordinates2 dup2=new Coordinates2(dupCover[cyc.get(a1+1)/2].dup.value,dupCover[cyc.get(a1+1)/2].dup.y,dupCover[cyc.get(a1+1)/2].dup.x,dupCover[cyc.get(a1+1)/2].dup.flag);	   
									    
									   Coordinates2 covdup2=new Coordinates2(dupCover[cyc.get(a1+1)/2].covDup.value,dupCover[cyc.get(a1+1)/2].covDup.y,dupCover[cyc.get(a1+1)/2].covDup.x,dupCover[cyc.get(a1+1)/2].covDup.flag);	   
										   
									   cov2.setCovDup(dup2,covdup2); // recuperer la deuxieme duplication
									   
									     
									   
								   }
								   
								   
								   
								   int seg=coorChev.value;
								   
								   //System.out.println("seg = "+seg);
								   
								   
								  // l'instruction en bas a ete remplacee par le bloc en haut
								   
								  //int seg=chev[v2/2][cyc.get(sauv+1)/2].value;
								   
								   
								   // fin de la partie ajoutee le 2 dec
								   
								   
								   
								   // mettre a jour la duplication et la partie qui la chevauche
								   
								   
								   dupCover[v2/2].updateDupCov(chev[v2/2][cyc.get(sauv+1)/2]);
								   
								   // Je dois augmenter le cout
								   
								   int reste=dupCover[v2/2].dup.value;
								   
								   //System.out.println("reste = "+reste);
								   
								   //System.out.println("reste "+reste+", segment = "+seg);
								   
								   
								   if (reste==0) 
								 
								   {
									   // debut de la partie ajoutee le 24 nov
									   
									   if (ctaille==5)
									   
									   {
										   
										   //System.out.println("Cycle de taille 5");
										   
										   //int a1=0;
										   
										   //while ( (a1<ctaille-1) && (cyc.get(a1)% 2 != 0) ) a1++;
										   
										   
										   if (cov3.dup.flag!=cov3.covDup.flag)  // source et cible dans dff�rents g�nomes
										   
										   {
											   
											   
											
											   //System.out.println("plus petite se trouve dans "+coorChev.flag);
											   
											   if (cov1.dup.x!=cov2.covDup.x || cov1.dup.y!=cov2.covDup.y) // ce n'est pas une transposition
											   
											   {
												   
												   /*
												   
												   System.out.println("cov1.dup.y "+cov1.dup.y);
												   
												   System.out.println("cov1.dup.x "+cov1.dup.x);
												   
												   System.out.println("cov1.covDup.y "+cov1.covDup.y);
												   
												   System.out.println("cov1.covDup.y "+cov1.covDup.x);
												   
												   System.out.println("cov2.dup.y "+cov2.dup.y);
												   
												   System.out.println("cov2.dup.x "+cov2.dup.x);
												   
												   System.out.println("cov2.covDup.y "+cov2.covDup.y);
												   
												   System.out.println("cov2.covDup.y "+cov2.covDup.x);
												   
												   
												   System.out.println("le pb est ici 1");
												   
												   
												   */
												   
												   if (flag==0)  // la plus petite se trouve dans X
												   
												   	{
													   
													   int type=isCoveredX[coorChev.y].dup.value;
													   
													   
													   for (int h=coorChev.y; h<=coorChev.x;h++)   
														   
													   {
														   
														   isCoveredX[h].dup.value=4;       // considerer chaque duplication comme perte dans Y
														   
														   isCoveredX[h].dup.y=h;          // position d'une perte entre (h,h)
														   isCoveredX[h].dup.x=h;
														   
														   // nbLY++; // incrementer le nombre de pertes
														   
														   
													   }
												   
													   /*
													   
													   System.out.println("type = "+type);
													   
													   
													   
													   System.out.println("cov1.dup.val = "+cov1.dup.value);
													   
													   
													   
													   System.out.println("cov1.covDup.val = "+cov1.covDup.value);
													   
													   
													   
													   System.out.println("cov2.dup.val = "+cov2.dup.value);
													   
													   
													   
													   System.out.println("cov2.covDup.val = "+cov2.covDup.value);
													   
													   
													   */
													   
													   nbLY=nbLY+seg;
													   
													   
													   
														   
													   if (type==7) 
														   
														   nbDupXY-- ;
													   
													   		else nbDupIXY-- ;
													  
													   
													   
												   	}
													   
												   
												   
												   else   // elle se trouve dans Y
												   
												   {
													   
													   int type=isCoveredY[coorChev.y].dup.value;
													   
													   for (int h=coorChev.y; h<=coorChev.x;h++)
														   
													   {
														   
														   isCoveredY[h].dup.value=9;     // considerer chaque duplication comme perte dans X
														   
														   isCoveredY[h].dup.y=h;        // position d'une perte entre (h,h)
														   isCoveredY[h].dup.x=h;
														   
														   // nbLX++; // incrementer le nombre de pertes
														   
													   }
												   
													   //System.out.println("Nombre de LX avant"+nbLX);
													   
												       nbLX=nbLX+seg;
												       
												       //System.out.println("Nombre de LX apres"+nbLX);
												       
												       
													   if (type==12) 
														   
														   nbDupYX-- ;
													   
													   		else nbDupIYX-- ;
													   
													   
												   }
												   
												   cout=cout-1+seg;   
												   
												   
											   }  
												   
											   
											   
											   else
												   
												   
											   {
												   
											   cout--;
											   nbTransp++;
											   
											   int type;
											   
											   if (flag==0) type=isCoveredX[coorChev.y].dup.value ;
											   else type=isCoveredY[coorChev.y].dup.value ;
											   
											   
											   
											   if (cov1.dup.flag==0)  // duplication XY ou duplication inverse XY
												   
											   {
												   
												   
												   for (int h=cov1.dup.y; h<=cov1.dup.x;h++)
													   
													   isCoveredX[h].dup.value=2;
											   
												   for (int h=cov2.dup.y; h<=cov2.dup.x;h++)
													   isCoveredY[h].dup.value=2;
												 
												   
											   }
											   
												   
											   else    
											   
											   {
												   
												   
												   for (int h=cov1.dup.y; h<=cov1.dup.x;h++)
													   isCoveredY[h].dup.value=2;
											   
												   
												   for (int h=cov2.dup.y; h<=cov2.dup.x;h++)
													   isCoveredX[h].dup.value=2;
												   
												   
												   
											   }
											   
											   if (flag==0) 
												   
											   	{
												   
												   if (type==7) 
													   
												   {
													   nbDupXY--;
													   nbDupYX--;
													   
												   }
												   
												   else 
													   
												   {
													   nbDupIXY--;
													   nbDupIYX--;
													   
												   }
													   
											   
											   
											   	}
											   
											   
											   else 
												   
											   	{
												   
												   if (type==12) 
													   
												   {
													   nbDupXY--;
													   nbDupYX--;
													   
												   }
												   
												   else 
													   
												   {
													   nbDupIXY--;
													   nbDupIYX--;
													   
												   }
													   
											   
											   
											   	}
											   
											   
											   
											   }
											   
											   
											   
										   }
										   
										   
										   else  // cycle de taille 5 dans le meme genome
											   
										   {
											  
											   
											   if (cov1.dup.flag==0)  // si la duplication est dans le genome X
											   
											   {
												   int typ=isCoveredX[coorChev.y].dup.value;
												   
												   
												   
												   
												   
												   for (int h=cov1.dup.y; h<=cov1.dup.x;h++)   
													   
												   {
													   
													   isCoveredX[h].dup.value=4;       // considerer chaque duplication comme perte dans Y
													   
													   isCoveredX[h].dup.y=h;          // position d'une perte entre (h,h)
													   isCoveredX[h].dup.x=h;
													   
													   // nbLY++; // incrementer le nombre de pertes
													   
													   
												   }
											   
												   
												   nbLY=nbLY+seg;
												   
												     
													   if (typ==5)
													   
												             nbDupX--;
												   
													   else 
														    nbDupIX--;
													   
													   
												   
											   }
											   
											   else    // sinon la duplication est dans le genome Y
											    
											   {
												   
												   int typ=isCoveredY[coorChev.y].dup.value;
												   
												   for (int h=cov1.dup.y; h<=cov1.dup.x;h++)
													   
												   {
													   
													   isCoveredY[h].dup.value=9;     // considerer chaque duplication comme perte dans X
													   
													   isCoveredY[h].dup.y=h;        // position d'une perte entre (h,h)
													   isCoveredY[h].dup.x=h;
													   
													   // nbLX++; // incrementer le nombre de pertes
													   
												   }
											   
												   //System.out.println("Nombre de LX avant"+nbLX);
												   
											       nbLX=nbLX+seg;
											       
											       //System.out.println("Nombre de LX apres"+nbLX);
											       
											       
											       
											       if (typ==10)
													   
											             nbDupY--;
											   
												   else 
													    nbDupIY--;
												   
												   
											   }
											   
											   
											   cout=cout-1+seg;
										   }
										   
										   
										   
										   
									   }
									
									   else   // la taille du cycle est superieure a 5 mais le reste est 0
									   
										   
									   {
										   
										   
										   
										   if (cov3.dup.flag==0)   // si la duplication est dans le genome X
											   
										   {
											   
											   int typ =isCoveredX[coorChev.y].dup.value;
											   
											   for (int h=cov3.dup.y; h<=cov3.dup.x;h++)   
												   
											   {
												   
												   isCoveredX[h].dup.value=4;       // considerer chaque duplication comme perte dans Y
												   
												   isCoveredX[h].dup.y=h;          // position d'une perte entre (h,h)
												   isCoveredX[h].dup.x=h;
												   
												   // nbLY++; // incrementer le nombre de pertes
												   
											   }
											   
											   
											   nbLY=nbLY+seg;
											   
											   if (typ==5) nbDupX--;
											   
											   else if (typ==6) nbDupIX--;
											   
											   else if (typ==7) nbDupXY--;
											   
											   else if (typ==8) nbDupIXY--;
										   
										   }
										   
										   
										   
										   else   // sinon la duplication est dans le genome Y
										   
										   {
											   
											   int typ =isCoveredY[coorChev.y].dup.value;
											   
											   for (int h=cov3.dup.y; h<=cov3.dup.x;h++)
												   
											   {
												   
												   isCoveredY[h].dup.value=9;     // considerer chaque duplication comme perte dans X
												   
												   isCoveredY[h].dup.y=h;        // position d'une perte entre (h,h)
												   isCoveredY[h].dup.x=h;
												   
												   //  nbLX++; // incrementer le nombre de pertes
												   
											   }  
										   
											   //System.out.println("Nombre de LX avant"+nbLX);
											   
											   nbLX=nbLX+seg;
											   
											   //System.out.println("Nombre de LX apres"+nbLX);
											   
											   
											   if (typ==10) nbDupY--;
											   
											   else if (typ==11) nbDupIY--;
											   
											   else if (typ==12) nbDupYX--;
											   
											   else if (typ==13) nbDupIYX--;
											   
											   
											   
											   
										   } 
										   
									   cout=cout-1+seg;
									   
									   
									   }
									   
									   
									   
									   
								   }
								   
								   
								   
								   else  // le reste n'est pas 0
									   
								   	{   
									   
									   
									   
									   
									   	// debut de la partie ajoutee le 3dec
									   
									   if (coorChev.flag==0)     // si la duplication est dans le genome X
									   
									   {
										   //int typ =isCoveredX[coorChev.y].dup.value;
										   
										   
										   for (int h=coorChev.y; h<=coorChev.x;h++)   
											   
										   {
											   
											   isCoveredX[h].dup.value=4;       // considerer chaque duplication comme perte dans Y
											   
											   isCoveredX[h].dup.y=h;          // position d'une perte entre (h,h)
											   isCoveredX[h].dup.x=h;
											   
											   // nbLY++; // incrementer le nombre de pertes
											   
										   }
									   
									        
										   for (int h=dupCover[v2/2].dup.y; h<=dupCover[v2/2].dup.x; h++)
										   
										   {
											   isCoveredX[h].dup.y=dupCover[v2/2].dup.y;   // modifier les coordonnees de la duplication
											   isCoveredX[h].dup.x=dupCover[v2/2].dup.x;
										   
										   }
										   
										   
										   nbLY=nbLY+seg;
										   
										   
										   //if (typ==5) nbDupX--;
										   
										   //else if (typ==6) nbDupIX--;
										   
										   //else if (typ==7) nbDupXY--;
										   
										   //else if (typ==8) nbDupIXY--;
										   
										   
										   
										   
										   
									   }
									   
									   
									   else   // sinon la duplication est dans le genome Y
										   
									   {
										   
										   //int typ =isCoveredY[coorChev.y].dup.value;
										   
										   
										   //System.out.println("taille du segment = "+seg);
										   
										   //System.out.println("debut = "+coorChev.y);
										   
										   //System.out.println("fin = "+coorChev.x);
										   
										   for (int h=coorChev.y; h<=coorChev.x;h++)   
											   
										   {
											   
											   isCoveredY[h].dup.value=9;       // considerer chaque duplication comme perte dans Y
											   
											   isCoveredY[h].dup.y=h;          // position d'une perte entre (h,h)
											   isCoveredY[h].dup.x=h;
											   
											   //nbLX++; // incrementer le nombre de pertes
											   
										   }
									   
										   
										   for (int h=dupCover[v2/2].dup.y; h<=dupCover[v2/2].dup.x; h++)
											   
										   {
											   isCoveredY[h].dup.y=dupCover[v2/2].dup.y;    // modifier les coordonnees de la duplication
											   isCoveredY[h].dup.x=dupCover[v2/2].dup.x;
										   
										   }
										   
										   //System.out.println("Nombre de LX avant"+nbLX);
										   nbLX=nbLX+seg;
										   //System.out.println("Nombre de LX apres"+nbLX);
									   
										   
										   //if (typ==10) nbDupY--;
										   
										   //else if (typ==11) nbDupIY--;
										   
										   //else if (typ==12) nbDupYX--;
										   
										   //else if (typ==13) nbDupIYX--;
										   
										   
										   
									   }
									   
										   
										   
									   
									   
									   // fin de la partie ajoutee le 3dec
									   
									   cout=cout+seg;
								   
								   
								   	}
								   
								   
								   
								   //System.out.println("taille du seg : "+seg+", Le reste : "+reste+ ", valeur cout : "+cout);
								   
								   
								   for (int v3=0; v3<v2/2; v3++)
								   		{
									   
									   chev[v2/2][v3]=Chevauch2(dupCover[v2/2].dup, dupCover[v3].covDup);
									   chev[v3][v2/2]=Chevauch2(dupCover[v3].dup, dupCover[v2/2].covDup);
									   
									   if (chev[v2/2][v3].getValue()==0) gam.removeEdge(v2, 2*v3);	
								   	   
									   if (chev[v3][v2/2].getValue()==0) gam.removeEdge(2*v3+1, v2-1);
								   		
								   		}
								   
								   
								   for (int v3=v2/2+1; v3<=var; v3++)
							   		{
								   
								   chev[v2/2][v3]=Chevauch2(dupCover[v2/2].dup, dupCover[v3].covDup);
								   chev[v3][v2/2]=Chevauch2(dupCover[v3].dup, dupCover[v2/2].covDup);
								   
								   if (chev[v2/2][v3].getValue()==0) gam.removeEdge(v2, 2*v3);	
							   	   
								   if (chev[v3][v2/2].getValue()==0) gam.removeEdge(2*v3+1, v2-1);
								   
							   		}
								   
								   
								   
								   
							  	}
						    		
					  	}					
					  						
					  
					  
				  	}
				  
				 
			
			}
			
			
			DirectedCycle dc=new DirectedCycle(gam);
			  
			
			
			  //if (dc.hasCycle()) System.out.println("il existe toujoutrs des cycles");
			  
			  //else System.out.println("Tous les cycles ont ete elimines");
			
			
			}
			
			
			
			
			
			cout_acyc=cout;
			
			
			
			
			time = System.currentTimeMillis() - time;
			
			
			
			// debut du post-traitement du voisinage
			
			
			
			
			// fin du post-traitement du voisinage
			
			
			/* debut du pretraitement (succession de plusieurs pertes)
			
			for (int h=0; h<isCoveredX.length; h++)
			      
				if( h+1 <isCoveredX.length)
					{
						if (isCoveredX[h].dup.value==4 && isCoveredX[h+1].dup.value==4)
							
						{
						
							nbLY--;
							cout--;
							
						}
						
						
						
					}
			
			
			
			
			for (int h=0; h<isCoveredY.length; h++)
			      
				if( h+1 <isCoveredY.length)
					{
						if (isCoveredY[h].dup.value==9 && isCoveredY[h+1].dup.value==9)
							
						{
							
							nbLX--;
							cout--;
							
						}
						
						
						
					}
			
			
			
			// fin du pretraitement (succession de plusieurs pertes)
			
			*/
			
			//System.out.println("===cout de l'alignement acyclique=== : "+cout);
			
			//System.out.println(" Temps de calcul : " + time + " milliseconds");
			
			String format = String.format("%%0%dd", 2);  
			  
			String seconds = String.format(format, (time/1000) % 60);  
			String minutes = String.format(format, ((time/1000) % 3600) / 60);  
			String hours = String.format(format, (time/1000) / 3600); 
			
			hhmmss =  hours + ":" + minutes + ":" + seconds;
			
			//System.out.println(" Temps de calcul : "+hhmmss);
			
			
			/*
		  	
			System.out.println("nbMatch = "+nbMatch);
			
			
			System.out.println("nbLY = "+nbLY+", nbDupX = "+nbDupX+", nbDupIX = "+nbDupIX+", nbDupXY = "+nbDupXY+", nbDupIXY = "+nbDupIXY);
			
			System.out.println("nbLX = "+nbLX+", nbDupY = "+nbDupY+", nbDupIY = "+nbDupIY+", nbDupYX = "+nbDupYX+", nbDupIYX = "+nbDupIYX);
			
			System.out.println("nbInvXY = "+nbInvXY+", nbTransp = "+nbTransp+", nbSub = "+nbSub);
			
			System.out.println();
			
			System.out.println("Interpr�tation des g�nes de X");
			
			System.out.println();
			
			
			int h=0;
			
			while (h<isCoveredX.length)
				
				
			{
				
				if (isCoveredX[h].dup.value!=0)
			
					{
					    System.out.println("Les genes "+h+"... "+isCoveredX[h].dup.x+ ": Type couverture = "+isCoveredX[h].dup.value+", taille = "+isCoveredX[h].dup.y+", "+isCoveredX[h].dup.x);
					
					   h=isCoveredX[h].dup.x+1;

						   
						   
				    	//System.out.println("gene "+h+", Type couverture = "+isCoveredX[h].dup.value+", taille = "+isCoveredX[h].dup.y+", "+isCoveredX[h].dup.x);
			
					
					}
					
				else h++;
			
			}
			
			
			
			System.out.println();
			
			System.out.println("Interpr�tation des g�nes de Y");
			
			System.out.println();
			
			
			
			
			h=0;
			
			while (h<isCoveredY.length)
				
			{
				
				if (isCoveredY[h].dup.value!=0)
			
					{
					    System.out.println("Les genes "+h+"... "+isCoveredY[h].dup.x+ ": Type couverture = "+isCoveredY[h].dup.value+", taille = "+isCoveredY[h].dup.y+", "+isCoveredY[h].dup.x);
					
					   h=isCoveredY[h].dup.x+1;

						   
						   
				    	//System.out.println("gene "+h+", Type couverture = "+isCoveredX[h].dup.value+", taille = "+isCoveredX[h].dup.y+", "+isCoveredX[h].dup.x);
			
					
					}
				
				
				else h++;
			
			}
			
			
			
			*/
			
			
			
			
			
		  	}
		  
		  
		  
		  
		  public void Affichage()
			
			{
				
				System.out.println(">===cout de l'alignement == :"+cout);
				
				
				System.out.println();
				
				System.out.println(">===cout Olivier == : ");
				
				System.out.println();
				
				System.out.println(">Temps de calcul : " + time + " milliseconds");
				
				System.out.println();
				
				System.out.println("nbMatch = "+nbMatch);
				
				
				System.out.println("nbLY = "+nbLY+", nbDupX = "+nbDupX+", nbDupIX = "+nbDupIX+", nbDupXY = "+nbDupXY+", nbDupIXY = "+nbDupIXY);
				
				System.out.println("nbLX = "+nbLX+", nbDupY = "+nbDupY+", nbDupIY = "+nbDupIY+", nbDupYX = "+nbDupYX+", nbDupIYX = "+nbDupIYX);
				
				System.out.println("nbInvXY = "+nbInvXY+", nbTransp = "+nbTransp+", nbSub = "+nbSub);
				
				System.out.println();
				
				System.out.println("Interpr�tation des g�nes de X");
				
				System.out.println();
				
				
				
				int h=0;
				
				while (h<isCoveredX.length)
					
					
				{
					
					if (isCoveredX[h].dup.value!=0)
				
						{
							System.out.println("Les genes "+st1[h]+"... "+st1[isCoveredX[h].dup.x]+ ", Position "+h+"... "+isCoveredX[h].dup.x+ ": �v�nement : "+isCoveredX[h].dup.value);
						
						   h=isCoveredX[h].dup.x+1;

							   
					    	//System.out.println("gene "+h+", Type couverture = "+isCoveredX[h].dup.value+", taille = "+isCoveredX[h].dup.y+", "+isCoveredX[h].dup.x);
				
						}
						
					else h++;
				
				}
				  
				
				System.out.println();
				
				System.out.println("Interpr�tation des g�nes de Y");
				
				System.out.println();
				
				
				
				
				h=0;
				
				while (h<isCoveredY.length)
					
				{
					
					if (isCoveredY[h].dup.value!=0)
				
						{
							System.out.println("Les genes "+st2[h]+"... "+st2[isCoveredY[h].dup.x]+", Position "+h+"... "+isCoveredY[h].dup.x+ ": �v�nement = "+isCoveredY[h].dup.value);
						
						   h=isCoveredY[h].dup.x+1;

							   
							   
					    	//System.out.println("gene "+h+", Type couverture = "+isCoveredX[h].dup.value+", taille = "+isCoveredX[h].dup.y+", "+isCoveredX[h].dup.x);
				
						
						}
					
					
					else h++;
				
				}
				
				
			}
		  
		  
		  
		  
		  public void Affichage2() // m�me model de output que le programme linear de patrick
			
			{
				
				System.out.println(">=== Total cost = "+cout);
				
				
				System.out.println();
				
				
				
				
				
				
			}
		  
		  
		
		  
		  
		  public void post_traitement_voisin(int []voisin)
		  
		  	{
			  
			  Align al2=new Align (X, voisin);
				
				al2.createAlign();
				
				int h=0, t=this.isCoveredX.length;
				
				while (h<t)
					
					{
					
					   
					
						if (this.isCoveredX[h].dup.value==5 || this.isCoveredX[h].dup.value==6 || this.isCoveredX[h].dup.value==7 || this.isCoveredX[h].dup.value==8)  // segment de duplication
							
							if (this.isCoveredX[h].dup.x-this.isCoveredX[h].dup.y >= 1) // car toute dup de taille 1 est -- > perte
							
							{   
							   
							   
								int y=this.isCoveredX[h].dup.y;
							
								int x=this.isCoveredX[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false;
								
								
								while (i<=x && !bool)   
								
								{
									
								if (al2.isCoveredX[i].dup.value!=0)     // tester si tous les genes sont interpretes comme match dans voisin
								       bool=true;
								
								       else i++;
								
								}
								
								
								
								
								if (!bool)  // le segment est une perte dans Y, il faut donc le transformer en succession de pertes
								
								{
									
									if (this.isCoveredX[h].dup.value==5)
										nbDupX--;
										
										else if (this.isCoveredX[h].dup.value==6)
										  nbDupIX--;
										
										else if (this.isCoveredX[h].dup.value==7)
											  nbDupXY--;
										
										
										else if (this.isCoveredX[h].dup.value==8)
											  nbDupIXY--;
										
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredX[k].dup.value=4;
									
									   
									   this.isCoveredX[k].dup.y=k;
									   this.isCoveredX[k].dup.x=k;
									   
									   nbLY++;    // incrementer le nombre de perte dans Y
									   
									   cout++;   // incrementer le cout de l'alignement
									   
									}
									
									cout--;   // il faut enlever 1 si la taille de la duplication est >1
									
									
									
									
								}
								
									
							}
							
						
						
						
						h=this.isCoveredX[h].dup.x+1;   // passer au prochain gene
					
						
					}
				
				
				
				al2= new Align (Y, voisin);   // faire la meme chose pour Y
				
				al2.createAlign();
				
				h=0; t=this.isCoveredY.length;
				
				while (h<t)
					
					{
					
					
					
						if (this.isCoveredY[h].dup.value==10 || this.isCoveredY[h].dup.value==11 || this.isCoveredY[h].dup.value==12 || this.isCoveredY[h].dup.value==13)
							
							if (this.isCoveredY[h].dup.x-this.isCoveredY[h].dup.y >= 1)
							
							{
							
								int y=this.isCoveredY[h].dup.y;
							
								int x=this.isCoveredY[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false;
								
								
								while (i<=x && !bool)
									
								
								{	
									
								if (al2.isCoveredX[i].dup.value!=0)
								       bool=true;
								
								       else i++;
								
								}
								
								
								
								if (!bool)  // le segment est une perte, il faut donc le transformer
								
								{
									
									//System.out.println("match entre "+y+", "+x+" de Y"+"position dans voisin "+al2.isCoveredX[h].covDup.y+", "+al2.isCoveredX[h].covDup.x);
									
									if (this.isCoveredY[h].dup.value==10)
									nbDupY--;
									
									else if (this.isCoveredY[h].dup.value==11)
									  nbDupIY--;
									
									else if (this.isCoveredY[h].dup.value==12)
										  nbDupYX--;
									
									
									else if (this.isCoveredY[h].dup.value==13)
										  nbDupIYX--;
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredY[k].dup.value=9;
									
									   
									   this.isCoveredY[k].dup.y=k;
									   this.isCoveredY[k].dup.x=k;
									   
									   nbLX++;
									   
									   cout++;
									   
									}
									
									cout--;
									
									
									
									
								}
								
									
							}
							
					
						h=this.isCoveredY[h].dup.x+1;
					
					}
				
				
				
		  	}
		  
		  
		  
		  
		  
		  public void post_traitement_voisin_perte(int []voisin)
		  
		  	{
			  
			  
			  // Ceci est un post_traitement permettant de transformer les pertes en duplication.
			  // Ceci pourrait �tre possible �tant donn� que le nombre d'alphabet des ARNt est tr�s petit par rapport 
			  // � la taille du g�nome.
			  
			  //En cas d'utilisation, il ne faut pas prendre en compte la procedure dup-->perte ou bien il faut tenir compte du tableau
			  // dup_petre
			  
			  
			  
			  Align al2=new Align (X, voisin);
				
				al2.createAlign();
				
				int h=0, t=this.isCoveredX.length;
				
				while (h<t)
					
					{
					
						if (this.isCoveredX[h].dup.value==4 || this.isCoveredX[h].dup.value==5 || this.isCoveredX[h].dup.value==6 || this.isCoveredX[h].dup.value==7 || this.isCoveredX[h].dup.value==8)  
							
							{   
							   
							   
								int y=this.isCoveredX[h].dup.y;
							
								int x=this.isCoveredX[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false; // par defaut, on consid�re que les genes sont pr�sents dans le voisin
								
								
								while (i<=x && !bool)   
								
								{
									
								if (al2.isCoveredX[i].dup.value!=0)     // tester si tous les genes sont interpretes comme match dans voisin
								       bool=true;
								
								       else i++;
								
								}
								
								
								
								if (bool && this.isCoveredX[h].dup.value==4)   // transformer la perte en duplication
									
									{
									
									
										for (int k=y; k<=x; k++)
											
											{
												this.isCoveredX[k].dup.value=5;  // la valeur 5 indique une transformation P-->Dup
												this.isCoveredX[k].dup.y=k;
												this.isCoveredX[k].dup.x=k;
											
												dup_perteX[k]=false;
											
											}
									
										nbLY--;
										
										nbDupX++;
										
									}
								
								
								
								if (!bool && this.isCoveredX[h].dup.value!=4)  // le segment est une perte dans Y, il faut donc le transformer en succession de pertes
								
								{
									
									// if (this.isCoveredX[h].dup.value==4)  ne rien faire
										
									// nbDupX--;
									
									
									if (this.isCoveredX[h].dup.value==5)
										nbDupX--;
										
										else if (this.isCoveredX[h].dup.value==6)
										  nbDupIX--;
										
										else if (this.isCoveredX[h].dup.value==7)
											  nbDupXY--;
										
										
										else if (this.isCoveredX[h].dup.value==8)
											  nbDupIXY--;
										
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredX[k].dup.value=4;
									
									   
									   this.isCoveredX[k].dup.y=k;
									   this.isCoveredX[k].dup.x=k;
									   
									   nbLY++;    // incrementer le nombre de perte dans Y
									   
									   cout++;   // incrementer le cout de l'alignement
									   
									}
									
									cout--;   // il faut enlever 1 si la taille de la duplication est >1
									
									
								}
								
								
							}
							
						
						
						h=this.isCoveredX[h].dup.x+1;   // passer au prochain gene
					
						
					}
				
				
				
				al2= new Align (Y, voisin);   // faire la meme chose pour Y
				
				al2.createAlign();
				
				h=0; t=this.isCoveredY.length;
				
				while (h<t)
					
					{
					
					
					
						if (this.isCoveredY[h].dup.value==9 || this.isCoveredY[h].dup.value==10 || this.isCoveredY[h].dup.value==11 || this.isCoveredY[h].dup.value==12 || this.isCoveredY[h].dup.value==13)
							
							
							{
							
							
								int y=this.isCoveredY[h].dup.y;
							
								int x=this.isCoveredY[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false;
								
								
								while (i<=x && !bool)
									
								
								{	
									
								if (al2.isCoveredX[i].dup.value!=0)
								       bool=true;
								
								       else i++;
								
								}
								
								
								
								
								if (bool && this.isCoveredY[h].dup.value==9)   // transformer la perte en duplication
									
								{
									for (int k=y; k<=x; k++)
										
										{
											this.isCoveredY[k].dup.value=10;  // la valeur -10 indique une transformation P-->Dup
											this.isCoveredY[k].dup.y=k;
											this.isCoveredY[k].dup.x=k;
										
											dup_perteY[k]=false;
										
										}
								
									
									nbLX--;
									
									nbDupY++;
									
									
								}
								
								
								if (!bool && this.isCoveredY[h].dup.value!=9)  // le segment est une perte, il faut donc le transformer
								
								{
									
									//System.out.println("match entre "+y+", "+x+" de Y"+"position dans voisin "+al2.isCoveredX[h].covDup.y+", "+al2.isCoveredX[h].covDup.x);
									
									// if (this.isCoveredY[h].dup.value==9)  ne rien faire
									
									// nbDupY--;
									
									
									
									if (this.isCoveredY[h].dup.value==10)
									nbDupY--;
									
									else if (this.isCoveredY[h].dup.value==11)
									  nbDupIY--;
									
									else if (this.isCoveredY[h].dup.value==12)
										  nbDupYX--;
									
									
									else if (this.isCoveredY[h].dup.value==13)
										  nbDupIYX--;
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredY[k].dup.value=9;
									
									   
									   this.isCoveredY[k].dup.y=k;
									   this.isCoveredY[k].dup.x=k;
									   
									   nbLX++;
									   
									   cout++;
									   
									}
									
									cout--;
									
									
									
									
								}
								
									
							}
							
					
						h=this.isCoveredY[h].dup.x+1;
					
					}
				
				
				
		  	}
		  
		  
		 
		  
		  
		  public void post_traitement_dupVperte()
		  
		  	{
			  
			  // transformer les duplications de taille 1 en pertes
			  
			  // necessite l'execussion du poste traitement de voisinage avant.
			  
			  // j'ai ajout� les tableau dup_perteX et dup_perteY pour pouvoir d�cider si la transformation est possible ou pas.
			  
			  
			  for (int h=0; h<isCoveredX.length; h++)
			      
					
				  			if (isCoveredX[h].dup.y==isCoveredX[h].dup.x && isCoveredX[h].dup.value>=5 && isCoveredX[h].dup.value<=8 && dup_perteX[h])
				  				
				  			
				  			{
				  
							if (isCoveredX[h].dup.value==5)
								
							{
							
								nbDupX--;
								
							}
						
							
							else if (isCoveredX[h].dup.value==6)
								
							{
							
								nbDupIX--;
								
							}
							
							
							
							else if (isCoveredX[h].dup.value==7)
								
							{
							
								nbDupXY--;
								
							}
							
							
							else if (isCoveredX[h].dup.value==8)
								
							{
							
								//System.out.println("isCoveredX[h].dup.value = "+isCoveredX[h].dup.value+" h ="+h);
								
								//System.out.println("isCoveredX[h].dup.x = "+isCoveredX[h].dup.value);
								
								
								
								//System.out.println("je suis la = "+h);
								nbDupIXY--;
								
							}
							
							
							isCoveredX[h].dup.value=4;
							nbLY++;
							
				  			}	
				
				
				
			  
			  
			  
			  
				for (int h=0; h<isCoveredY.length; h++)
				
					
					
					if (isCoveredY[h].dup.y==isCoveredY[h].dup.x  && isCoveredY[h].dup.value>=10 && isCoveredY[h].dup.value<=13 && dup_perteY[h])
					
					{
						
						
					if (isCoveredY[h].dup.value==10)
						
					{
					
						nbDupY--;
						
					}
					
					
					else if (isCoveredY[h].dup.value==11)
						
					{
					
						nbDupIY--;
						
					}
					
					
					else if (isCoveredY[h].dup.value==12)
						
					{
					
						nbDupYX--;
						
					}
					
					

					else if (isCoveredY[h].dup.value==13)
						
					{
					
						nbDupIYX--;
						
					}
					
					
					isCoveredY[h].dup.value=9;
					
					nbLX++;
					
					
					}
					
					
						
		  	}
		  
		  
		  
		  
		  public void post_traitement_pertes_succ()
		  
		  	{
			  
			  // regrouper les pertes successives en une seule perte
			  
			  
			// debut du pretraitement (succession de plusieurs pertes)
				
				for (int h=0; h<isCoveredX.length; h++)
				      
					if( h+1 <isCoveredX.length)
						{
							if (isCoveredX[h].dup.value==4 && isCoveredX[h+1].dup.value==4)
								
							{
							
								nbLY--;
								cout--;
								
							}
							
							
							
						}
				
				
				
				
				for (int h=0; h<isCoveredY.length; h++)
				      
					if( h+1 <isCoveredY.length)
						{
							if (isCoveredY[h].dup.value==9 && isCoveredY[h+1].dup.value==9)
								
							{
								
								nbLX--;
								cout--;
								
							}
							
							
							
						}
				
				
				
				// fin du pretraitement (succession de plusieurs pertes)
			  
			  
			  
		  	}
		  
		  
		  
		  
		  public void pertes_succ()
		  
		  	{
			  
			  // mettre � jour les taills des segments de pertes apr�s leurs regroupements
			  
			  int h=0;
			
			  while (h<isCoveredX.length)
			  	
			  	{
				  
				  if (isCoveredX[h].dup.value==4)
					  
				  	{
					  
					  int debut=h;
					  	
					  while (debut<isCoveredX.length && isCoveredX[debut].dup.value==4)
						  	debut++;
					  
					 isCoveredX[h].dup.x=debut-1;
					  
					  
				  	}
				  
				  h=isCoveredX[h].dup.x+1;
				  
			  	}
			  
			  
			  
			  h=0;
				
			  while (h<isCoveredY.length)
			  	
			  	{
				  
				  if (isCoveredY[h].dup.value==9)
					  
				  	{
					  
					  int debut=h;
					  	
					  while (debut<isCoveredY.length && isCoveredY[debut].dup.value==9)
						  	debut++;
					  
					  isCoveredY[h].dup.x=debut-1;
					  
					  
				  	}
				  
				  h=isCoveredY[h].dup.x+1;
				  
			  	}
			  
				
				
				// fin du pretraitement (succession de plusieurs pertes)
			  
			  
			  
		  	}
		  
		  
		  
		  
		  public void post_traitement_voisin2(int []voisin) // prendre en compte les substitutions, transpositions et inversions
		  
		  	{
			  
			  Align al2=new Align (X, voisin);   // commencer par X
				
				al2.createAlign();
				
				int h=0, t=this.isCoveredX.length;
				
				while (h<t)
					
					{
					
					   
					
						if (this.isCoveredX[h].dup.value==1 || this.isCoveredX[h].dup.value==2 || this.isCoveredX[h].dup.value==3 || this.isCoveredX[h].dup.value==5 || this.isCoveredX[h].dup.value==6 || this.isCoveredX[h].dup.value==7 || this.isCoveredX[h].dup.value==8)  // segment de duplication
							
							
							{   
							   
							   
								int y=this.isCoveredX[h].dup.y;
							
								int x=this.isCoveredX[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false;
								
								
								while (i<=x && !bool)   
								
								{
									
								if (al2.isCoveredX[h].dup.value!=0)     // tester si tous les genes sont interpretes comme match dans voisin
								       bool=true;
								
								       else i++;
								
								}
								
								
								
								
								if (!bool)  // le segment est une perte dans Y, il faut donc le transformer en succession de pertes 
								           //  ou bien choisir l'inversion ou la substitution dans l'autre genome
								{
								
									if (this.isCoveredX[h].dup.value==1) //substitution
										
										for (int k=y; k<=x; k++)
										   this.isCoveredX[k].dup.value=14;
									
									
									else if (this.isCoveredX[h].dup.value==2) //transposition
										
										for (int k=y; k<=x; k++)
											
										   this.isCoveredX[k].dup.value=15;
									
										
									else if (this.isCoveredX[h].dup.value==3) //inversion
										
										for (int k=y; k<=x; k++)
											
										   this.isCoveredX[k].dup.value=16;
									
									
									
									
									
								else
									
								{
									
									
									if (this.isCoveredX[h].dup.value==5)
										nbDupX--;
										
										else if (this.isCoveredX[h].dup.value==6)
										  nbDupIX--;
										
										else if (this.isCoveredX[h].dup.value==7)
											  nbDupXY--;
										
										
										else if (this.isCoveredX[h].dup.value==8)
											  nbDupIXY--;
										
									
									int taille=this.isCoveredX[h].dup.x-this.isCoveredX[h].dup.y+1;
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredX[k].dup.value=4;
									
									   
									   this.isCoveredX[k].dup.y=k;
									   this.isCoveredX[k].dup.x=k;
									   
									   nbLY++;    // incrementer le nombre de perte dans Y
									   
									   cout++;   // incrementer le cout de l'alignement
									   
									}
									
									if (taille >1 ) cout--;   // il faut toujours enlever 1 du cout si la duplication est de taille > 1
									
									
								}	
								
								
								
								}
								
								
								
								
								
								
								
								
							}
							
						
						
						
						h=this.isCoveredX[h].dup.x+1;   // passer au prochain gene
					
						
					}
				
				
				
				al2= new Align (Y, voisin);   // faire la meme chose pour Y
				
				al2.createAlign();
				
				h=0; t=this.isCoveredY.length;
				
				while (h<t)
					
					{
					
					
					
						if (this.isCoveredY[h].dup.value==10 || this.isCoveredY[h].dup.value==11 || this.isCoveredY[h].dup.value==12 || this.isCoveredY[h].dup.value==13)
							
							
							{
							
								int y=this.isCoveredY[h].dup.y;
							
								int x=this.isCoveredY[h].dup.x;
								
								
								int i=y;
								
								boolean bool=false;
								
								
								while (i<=x && !bool)
								
								if (al2.isCoveredX[h].dup.value!=0)
								       bool=true;
								
								       else i++;
								
								
								if (!bool)  // le segment est une perte, il faut donc le transformer
								
								{
									
									
									
									if (this.isCoveredY[h].dup.value==10)
									nbDupY--;
									
									else if (this.isCoveredY[h].dup.value==11)
									  nbDupIY--;
									
									else if (this.isCoveredY[h].dup.value==12)
										  nbDupYX--;
									
									
									else if (this.isCoveredY[h].dup.value==13)
										  nbDupIYX--;
									
									
									int taille=this.isCoveredY[h].dup.x-this.isCoveredY[h].dup.y+1;
									
									
									for (int k=y; k<=x; k++)
										
									{
										
									   this.isCoveredY[k].dup.value=9;
									
									   
									   this.isCoveredY[k].dup.y=k;
									   this.isCoveredY[k].dup.x=k;
									   
									   nbLX++;
									   
									   cout++;
									   
									}
									
									if (taille >1 ) cout--; // decrementer si la taille de la duplication >1
									
									
									
									
								}
								
									
							}
							
					
						h=this.isCoveredY[h].dup.x+1;
					
					}
				
				
				
				
				
				
				
				
				
				
		  	}
		  
		  
		  
		  public void AlignOut1()   // Afficher l'alignement des deux s�quences
		  
		  
		  	{
			  
			  
			  
			  String space="   ---   ";
			  
			  String vide="   ";
			  
			  String point=".";
			  
			  String formatX="";
			  
			  String formatY="";
			  
			  int Xlength=isCoveredX.length, Ylength=isCoveredY.length;
				
				int i=0,j=0;
			 
			  String []st3=this.st1.clone();
				
			  String []st4=this.st2.clone();
			  
			  
			  formatX=String.format("%11s","Genome X");
		   		
		   	  formatY=String.format("%-11s","Genome Y");
		   
		   	System.out.println();
		   	  System.out.println("Detailed Results");
		   	 System.out.println("****************");
			  
			  System.out.println();
		   	  
		   	  System.out.println(formatX+space+formatY+vide+"Operation");
			  
		   	  formatX=String.format("%11s","========");
	   		
		   	  formatY=String.format("%-11s","========");
		   
		   	  System.out.println(formatX+space+formatY+vide+"=========");
		   	  
		   	  System.out.println();
			  
			  
			 
			  for (int l=0; l<Xlength; l++)
				  
			  
				  st3[l]=st1[l]+"_"+(l+1);
			  
			  
			  for (int l=0; l<Ylength; l++)
			
				  st4[l]=st2[l]+"_"+(Xlength+l+1);
			   
				  
				
				while (i<Xlength && j<Ylength)
					
					
					{
						if (isCoveredX[i].dup.value<=3)  // match, substitution, transposition ou inversion
					       
							{
							
							   if (isCoveredX[i].dup.value==2) // transposition
							
							   		{
								   
								       
								   
								       formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   
								       formatY=String.format("%-11s",point);
								       
								       
								   
								   	   System.out.println(formatX+space+formatY+vide+"transposition of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" with "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
								   	  
								   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
							   			
						   				{
								   		   
								   		formatX=String.format("%11s",st3[h]);
										   
								   		System.out.println(formatX+space+formatY);
								   		   
						   				}
								   	   
								   	   
								   	   
								       i=isCoveredX[i].dup.x+1;
								       
							   		}
								   
							   
							   else   // match ou inversion ou substitution
								   
								   
							   	{
								   
								   while (isCoveredY[j].dup.value!=0 && isCoveredY[j].dup.value!=1 && isCoveredY[j].dup.value!=3)
								   
								   		{
									   
									       if (isCoveredY[j].dup.value==9)  // LX
									    	   
									       {
											   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"Loss of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	   
										   	   
									   		}
									      
									   
									       else if (isCoveredY[j].dup.value==10)  // Dy
									    	   
									       {
											   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st4[isCoveredY[j].covDup.y]+"..."+st4[isCoveredY[j].covDup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	   
									   		}
									       
									       
									       
									       else if (isCoveredY[j].dup.value==11)  // DIy
									    	   
									       {
											   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st4[isCoveredY[j].covDup.y]+"..."+st4[isCoveredY[j].covDup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	   
										   	   
									   		}
									       
									       
									       else if (isCoveredY[j].dup.value==12)  // Dyx
									    	   
									       {
											   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	   
										   	   
									   		}
									       
									       
									       
									       else if (isCoveredY[j].dup.value==13)  // DIyx
									    	   
									       {
											   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	   
										       
									   		}
									       
									       
									       else if (isCoveredY[j].dup.value==2) // transposition
												
									   		{
										   
										       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
										   
										       formatX=String.format("%11s",point);
										       
										       
										   
										   	   System.out.println(formatX+space+formatY+vide+"transposition of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" with "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
										   	  
										   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
									   			
								   				{
										   		   
										   		formatY=String.format("%-11s",st4[h]);
												   
										   		System.out.println(formatX+space+formatY);
										   		   
								   				}
										   	   
										   	  
									   		}
									       
									   
									       j=isCoveredY[j].dup.x+1;
									   		
									   			
								   		}
									
								   
								   
								   if (isCoveredX[i].dup.value==1) // substitution
									   
								   		{
									   
									    
									   
									    
									   		formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
									   		
									   		formatY=String.format("%-11s",st4[isCoveredX[i].covDup.y]);
									   
									   		System.out.println(formatX+space+formatY+vide+"Substitution");
									   		
									   
									   		
									   		
								   		}
								   
								   
								   else if (isCoveredX[i].dup.value==3) // inversion
									   
							   		{
								   
								   		formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   		
								   		formatY=String.format("%-11s",st4[isCoveredX[i].covDup.y]);
								   
								   		System.out.println(formatX+space+formatY+vide+"Inversion of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" with "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);
								   		
								        int iv=isCoveredX[i].covDup.y+1; // pour afficher les genes invers�s de Y
								   		
								   		
								   		for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
								   			
						   				{
								   		   
								   			formatX=String.format("%11s",st3[h]);
									   		
									   		formatY=String.format("%-11s",st4[iv]);
									   		
										   
								   		System.out.println(formatX+space+formatY);
								   		   
								   			iv++;
								   		
						   				}
								   		
								   		
								   		
							   		}
							   
								   
								   
								   else // commencer le match entre X et Y
								   
								   {
								   
									   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   		
								   	   formatY=String.format("%-11s",st4[isCoveredX[i].covDup.y]);
								   
								   	   System.out.println(formatX+space+formatY+vide);
									   
								   		
								   
								   }
								   	
								   
								   
								   
								   i=isCoveredX[i].dup.x+1;
							   	   j=isCoveredY[j].dup.x+1;
								   
							   	   
								   		
							   	}
							
							   
							
							}
						
						
						
						else
							
							{
							
							// commencer avec X
							
							
							if (isCoveredX[i].dup.value==4)  // LY
						    	   
						       {
								   
							       formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
							   
							       formatY=String.format("%-11s",point);
							       
							      
							   
							   	   System.out.println(formatX+space+formatY+vide+"Loss of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]);	
							   	  
							   	   
							   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
						   			
					   				{
							   		   
							   		formatX=String.format("%11s",st3[h]);
									   
							   		System.out.println(formatX+space+formatY);
							   		   
					   				}
							   	   
							   	   
						   		}
						      
						   
						       else if (isCoveredX[i].dup.value==5)  // Dx
						    	   
						       {
								   
						    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   
							       formatY=String.format("%-11s",point);
							       
							       
							   
							   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st3[isCoveredX[i].covDup.y]+"..."+st3[isCoveredX[i].covDup.x]);	
							   	  
							   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
						   			
					   				{
							   		   
							   		formatX=String.format("%11s",st3[h]);
									   
							   		System.out.println(formatX+space+formatY);
							   		   
					   				}
							   	   
							   	  
						   		}
						       
						       
						       
						       else if (isCoveredX[i].dup.value==6)  // DIx
						    	   
						       {
								   
						    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   
							       formatY=String.format("%-11s",point);
							       
							       
							   
							   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st3[isCoveredX[i].covDup.y]+"..."+st3[isCoveredX[i].covDup.x]);	
							   	  
							   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
						   			
					   				{
							   		   
							   		formatX=String.format("%11s",st3[h]);
									   
							   		System.out.println(formatX+space+formatY);
							   		   
					   				}
							   	   
							   	   
						   		}
						       
						       
						       else if (isCoveredX[i].dup.value==7)  // Dxy
						    	   
						       {
								   
						    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   
							       formatY=String.format("%-11s",point);
							       
							       
							   
							   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
							   	  
							   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
						   			
					   				{
							   		   
							   		formatX=String.format("%11s",st3[h]);
									   
							   		System.out.println(formatX+space+formatY);
							   		   
					   				}
							   	   
							   	   
							   	   
						   		}
						       
						       
						       
						       else if (isCoveredX[i].dup.value==8)  // DIxy
						    	   
						       {
								   
						    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
								   
							       formatY=String.format("%-11s",point);
							       
							       
							   
							   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
							   	  
							   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
						   			
					   				{
							   		   
							   		formatX=String.format("%11s",st3[h]);
									   
							   		System.out.println(formatX+space+formatY);
							   		   
					   				}
							   	   
							   	   
						   		}
							
							
						       else if (isCoveredX[i].dup.value==2) // transposition
								
					   		{
						   
						       formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
						   
						       formatY=String.format("%-11s",point);
						       
						       
						   
						   	   System.out.println(formatX+space+formatY+vide+"transposition of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" with "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
						   	  
						   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
					   			
				   				{
						   		   
						   		formatX=String.format("%11s",st3[h]);
								   
						   		System.out.println(formatX+space+formatY);
						   		   
				   				}
						   	  
						       
					   		}
								 
								 
							 i=isCoveredX[i].dup.x+1;
							
							 
							}
						
						
					}
				
				
				//if (i<Xlength) System.out.println("Erreur");
				
				
				while (i<Xlength)
					
					
				{
					
					
					// commencer avec X
					
					
					if (isCoveredX[i].dup.value==4)  // LY
				    	   
				       {
						   
					       formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
					   
					       formatY=String.format("%-11s",point);
					       
					      
					   
					   	   System.out.println(formatX+space+formatY+vide+"Loss of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]);	
					   	  
					   	   
					   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
				   			
			   				{
					   		   
					   		formatX=String.format("%11s",st3[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
				   		}
				      
				   
				       else if (isCoveredX[i].dup.value==5)  // Dx
				    	   
				       {
						   
				    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
						   
					       formatY=String.format("%-11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st3[isCoveredX[i].covDup.y]+"..."+st3[isCoveredX[i].covDup.x]);	
					   	  
					   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
				   			
			   				{
					   		   
					   		formatX=String.format("%11s",st3[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	  
				   		}
				       
				       
				       
				       else if (isCoveredX[i].dup.value==6)  // DIx      modifi� le 9 �tait isCoveredY, modif le 13 �tait [j]
				    	   
				       {
						   
				    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
						   
					       formatY=String.format("%-11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st3[isCoveredX[i].covDup.y]+"..."+st3[isCoveredX[i].covDup.x]);	
					   	  
					   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
				   			
			   				{
					   		   
					   		formatX=String.format("%11s",st3[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
				   		}
				       
				       
				       else if (isCoveredX[i].dup.value==7)  // Dxy   modifi� le 9 �tait isCoveredY, modif le 13 �tait [j]
				    	   
				       {
						   
				    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
						   
					       formatY=String.format("%-11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
					   	  
					   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
				   			
			   				{
					   		   
					   		formatX=String.format("%11s",st3[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
					   	   
				   		}
				       
				       
				       
				       else if (isCoveredX[i].dup.value==8)  // DIXy  // modif le 13 �tait isCoveredY[j]
				    	   
				       {
						   
				    	   formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
						   
					       formatY=String.format("%-11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" copied from "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
					   	  
					   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
				   			
			   				{
					   		   
					   		formatX=String.format("%11s",st3[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
				   		}
						 
						 
				       else if (isCoveredX[i].dup.value==2) // transposition
						
			   		{
				   
				       formatX=String.format("%11s",st3[isCoveredX[i].dup.y]);
				   
				       formatY=String.format("%-11s",point);
				       
				       
				   
				   	   System.out.println(formatX+space+formatY+vide+"transposition of "+st3[isCoveredX[i].dup.y]+"..."+st3[isCoveredX[i].dup.x]+" with "+st4[isCoveredX[i].covDup.y]+"..."+st4[isCoveredX[i].covDup.x]);	
				   	  
				   	   for (int h=isCoveredX[i].dup.y+1; h<=isCoveredX[i].dup.x; h++)
			   			
		   				{
				   		   
				   		formatX=String.format("%11s",st3[h]);
						   
				   		System.out.println(formatX+space+formatY);
				   		   
		   				}
				   	   
				   	   
			   		}
					
					
					
					 i=isCoveredX[i].dup.x+1;
					
						
					
					
				}
			
			
				
				while (j<Ylength)
				
					{
					
					 if (isCoveredY[j].dup.value==9)  // LX
				    	   
				       {
						   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Loss of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
					   	   
				   		}
				      
				   
				       else if (isCoveredY[j].dup.value==10)  // Dy
				    	   
				       {
						   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st4[isCoveredY[j].covDup.y]+"..."+st4[isCoveredY[j].covDup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
				   		}
				       
				       
				       
				       else if (isCoveredY[j].dup.value==11)  // DIy
				    	   
				       {
						   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st4[isCoveredY[j].covDup.y]+"..."+st4[isCoveredY[j].covDup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
					   	   
				   		}
				       
				       
				       else if (isCoveredY[j].dup.value==12)  // Dyx
				    	   
				       {
						   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
					   	   
				   		}
				       
				       
				       
				       else if (isCoveredY[j].dup.value==13)  // DIyx
				    	   
				       {
						   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"Inverted duplication of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" copied from "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	   
					       
				   		}
					 
					 
					 
				       else if (isCoveredY[j].dup.value==2) // transposition
							
				   		{
					   
					       formatY=String.format("%-11s",st4[isCoveredY[j].dup.y]);
					   
					       formatX=String.format("%11s",point);
					       
					       
					   
					   	   System.out.println(formatX+space+formatY+vide+"transposition of "+st4[isCoveredY[j].dup.y]+"..."+st4[isCoveredY[j].dup.x]+" with "+st3[isCoveredY[j].covDup.y]+"..."+st3[isCoveredY[j].covDup.x]);	
					   	  
					   	   for (int h=isCoveredY[j].dup.y+1; h<=isCoveredY[j].dup.x; h++)
				   			
			   				{
					   		   
					   		formatY=String.format("%-11s",st4[h]);
							   
					   		System.out.println(formatX+space+formatY);
					   		   
			   				}
					   	   
					   	  
				   		}
					 
				       
				   
				       j=isCoveredY[j].dup.x+1;
					
					}
			  
				
				System.out.println();
				System.out.println();
				
			  
		  	}
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  public void get_distribution_events()
		  
		  	{
			  
			  // determiner la distribution des �v�nements ( le nombre de chaque type d'�v�nements et sa taille )
			  
			  Arrays.fill(stat_dup, 0);
			  Arrays.fill(stat_loss, 0);
			  Arrays.fill(stat_inver, 0);
			  Arrays.fill(stat_trans, 0);
			  
			  
			  int h=0;
			
			  while (h<isCoveredX.length)
			  	
			  	{
				  
				  if (isCoveredX[h].dup.value==2) 
					  
				  {
					  
					  stat_trans[isCoveredX[h].dup.x-isCoveredX[h].dup.y+1]++;
					  
				  }
				  
				  
				  if (isCoveredX[h].dup.value==3)
				  {
						
					  stat_inver[isCoveredX[h].dup.x-isCoveredX[h].dup.y+1]++;
				  
					  
				  }
				  
				  
				  if (isCoveredX[h].dup.value==4) 
					  
				  {
					  stat_loss[isCoveredX[h].dup.x-isCoveredX[h].dup.y+1]++;
				  
					  
				  }
					  
				  if (isCoveredX[h].dup.value==5 || isCoveredX[h].dup.value==6 || isCoveredX[h].dup.value==7 || isCoveredX[h].dup.value==8)
					  
				  {
						  
					  stat_dup[isCoveredX[h].dup.x-isCoveredX[h].dup.y+1]++;
				  
				  }
				  
				  
				  h=isCoveredX[h].dup.x+1;
				  
			  	}
			  
			  
			  
			  h=0;
				
			  while (h<isCoveredY.length)
			  	
			  	{
				  // ne pas compter les transpositions et les inversions
				  
				  if (isCoveredY[h].dup.value==9) stat_loss[isCoveredY[h].dup.x-isCoveredY[h].dup.y+1]++;
				  
				  if (isCoveredY[h].dup.value==10 || isCoveredY[h].dup.value==11 || isCoveredY[h].dup.value==12 || isCoveredY[h].dup.value==13)
					  stat_dup[isCoveredY[h].dup.x-isCoveredY[h].dup.y+1]++;
				  
				  h=isCoveredY[h].dup.x+1;
				  
			  	}
			  
				
				
				// fin du pretraitement (succession de plusieurs pertes)
			  
			  
			  
		  	}
		  
		  
		  
		  public int[][] getTablePD()
		  {
			  return this.tablePD;
		  }
		  
		  
		  
		  
		  public void printTablePD()
		  
		  	{
			  for (int j=0; j<Y.length+1; j++)
			  
			  {
				  for (int i=0; i<X.length+1; i++)
					  System.out.print(tablePD[j][i]+"  ");
				  
				  System.out.println();
			  
				  
			  }
			  
		  	}
		  
		  
		   long get_timeValue()
		   	{
			   return time;
			   
		   	}
		 
		  
		   String get_time()
		   	{
			   return hhmmss;
			   
		   	}
		  
		  
		   
		   int get_coutCyc()
		   	{
			   return cout_cyc;
		   	}
		   
		  
		   
		   int get_coutAcyc()
		   	{
			   return cout_acyc;
		   	}
		  
		   
		   boolean hasCycle()
		   		{
			   		return hascycle;		
			   					
		   		}
		   
		   
		   
		   /**
			 *
			 *  cDup: calculer le cout d'une duplication
			 */
		   
		   int cDup(int i)
		   {
			   return 1;
			   
		   }
		   
		   
		   
		 
		   /**
			 *
			 *  cLost: calculer le cout d'une perte
			 */  
		   
		 int cLost(int seg)
		 
		 {
			return seg; 
			 
			 
		 }
		   
		   
		 
		
	}