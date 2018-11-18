package fem;


import math.BandMat;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpBlockMat;
import math.SpMat;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class MechMatrix {

	private Calculator calc;

	public MechMatrix(){}

	public MechMatrix(Model model){

		this.calc=new Calculator(model);

	}

	public void setTmpMat(Model model){
		setTmpMat(model,false);
	}

	public void setTmpMat(Model model,boolean massNeeded){

		System.out.println(" Structural analysis... ");

		System.out.println();
		System.out.println(" Number of unknown non-fixed points : "+model.numberOfUnknownT);
		System.out.println();

		System.out.println(" Calculating stiffness matrix ...");
		Mat Se=new Mat(model.nElVert,model.nElVert);
		Mat Me=new Mat(model.nElVert,model.nElVert);
		
		int m,row,column,rowNodeNumber,colNodeNumber,ext=4;
		int[] nz=new int[model.numberOfUnknownT];
		SpMat Ss=new SpMat(model.numberOfUnknownT);
		SpMat Ms=new SpMat(model.numberOfUnknownT);
		
		for(int i=0;i<model.numberOfUnknownT;i++){
			Ss.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		
		if(massNeeded)
			Ms.row[i]=new SpVect(model.numberOfUnknownT,model.nNodNod);
		}

		
		model.bT=new Vect(model.numberOfUnknownT);

		for(int i=1;i<=model.numberOfElements;i++){

			if(!model.element[i].isThermal()) {continue;}

			int[] vertNumb=model.element[i].getVertNumb();


			if(massNeeded){

				Se=this.calc.Se(model,i,Me);


			}
			else{
				Se=this.calc.Se(model,i);

			}

			for(int j=0;j<model.nElVert;j++){

				rowNodeNumber=vertNumb[j];

				if(model.node[rowNodeNumber].is_T_known()) continue;

				row=model.T_unknownIndex[rowNodeNumber]-1;

				for(int k=0;k<model.nElVert;k++){

					colNodeNumber=vertNumb[k];


					//=========================

					if(model.node[colNodeNumber].is_T_known() ){
						double T;
					//	if(model.node[colNodeNumber].getMap()==0) 
							T=model.node[colNodeNumber].T;									
						/*else {
							T=model.node[model.node[colNodeNumber].getMap()].T;
						}
*/
										model.bT.el[row]-=Se.el[j][k]*T;
							
										continue;
					}

				

					column=model.T_unknownIndex[colNodeNumber]-1;


					if(column>row) continue;

					m=util.search(Ss.row[row].index,nz[row]-1,column);

					if(m<0)
					{	
						
						Ss.row[row].index[nz[row]]=column;															
						Ss.row[row].el[nz[row]]=Se.el[j][k];
						if(massNeeded){
							Ms.row[row].index[nz[row]]=column;	
							Ms.row[row].el[nz[row]]=Me.el[j][k];
						}

						nz[row]++;

						if(nz[row]==Ss.row[row].nzLength-1){
							Ss.row[row].extend(ext);
							if(massNeeded)
								Ms.row[row].extend(ext);
						}

					}

					else{

						Ss.row[row].addToNz(Se.el[j][k],m);
						if(massNeeded)
							Ms.row[row].addToNz(Me.el[j][k],m);


					}
				}

			}
		
		}
	
		


	Ss.sortAndTrim(nz);


		model.Ss=Ss.deepCopy();
		if(massNeeded){

			Ms.sortAndTrim(nz);	
			model.Ms=Ms.deepCopy();
		}



	}

	public void setStiffMat(Model model){
		setStiffMat(model,false);
	}

	public void setStiffMat(Model model,boolean massNeeded){


		System.out.println(" Structural analysis... ");

		System.out.println();
		System.out.println(" Number of unknown non-fixed points : "+model.numberOfUnknownU);
		System.out.println();

		System.out.println(" Calculating stiffness matrix ...");
		Mat[][] Ke=new Mat[model.nElVert][model.nElVert];
		Mat[][] Me=new Mat[model.nElVert][model.nElVert];
		Mat RiT=new Mat(),Rj;
		int m,row,column,rowNodeNumber,colNodeNumber,rowb=0,ext=4;
		int[] nz=new int[model.numberOfUnknownU];
		SpBlockMat BKs=new SpBlockMat(model.numberOfUnknownU,model.numberOfUnknownU,model.nNodNod);
		SpBlockMat MKs=new SpBlockMat();
		if(massNeeded){
			MKs=new SpBlockMat(model.numberOfUnknownU,model.numberOfUnknownU,model.nNodNod);
		}

		int[][] unknown_U_ind=new int[model.numberOfUnknownU][model.dim];
		int ix=1;
		for(int i=0;i<unknown_U_ind.length;i++){
			int nodeNumb=model.unknownUnumber[i+1];

			for(int k=0;k<model.dim;k++){			
				if(!model.node[nodeNumb].is_U_known(k)){
					unknown_U_ind[i][k]=ix++;
				}
			}
		}


		model.bU=new Vect(model.numberOfUnknownUcomp);

		for(int i=1;i<=model.numberOfElements;i++){

			if(!model.element[i].isDeformable()) {continue;}
			

			int[] vertNumb=model.element[i].getVertNumb();


			if(massNeeded){

				Ke=this.calc.Ke(model,i,Me);

			}
			else{
				Ke=this.calc.Ke(model,i);

			}

			for(int j=0;j<model.nElVert;j++){

				rowNodeNumber=vertNumb[j];

				if(model.node[rowNodeNumber].is_U_known()) continue;

				if(model.coordCode==1)
					if(model.dim==2){
						RiT=util.rotMat2D(util.getAng(model.node[rowNodeNumber].getCoord())).transp();

					}
					else {
						double alpha=util.getAng(model.node[rowNodeNumber].getCoord().v2());
						RiT=util.rotEuler(new Vect(0,0,1), alpha).transp();
						
					}

				row=model.U_unknownIndex[rowNodeNumber]-1;

				for(int k=0;k<model.nElVert;k++){

					colNodeNumber=vertNumb[k];

					if(model.coordCode==1){
					

						if(model.dim==2){
							Rj=util.rotMat2D(util.getAng(model.node[colNodeNumber].getCoord()));
						}
						else {
							
							double alpha=util.getAng(model.node[colNodeNumber].getCoord().v2());
							Rj=util.rotEuler(new Vect(0,0,1), alpha);
						
						
						}

					
						Ke[j][k]=RiT.mul(Ke[j][k].mul(Rj));
						
						if(massNeeded)
							Me[j][k]=RiT.mul(Me[j][k].mul(Rj));
					}

					//=========================

					if(model.node[colNodeNumber].has_U_known() ){
						Vect v;
						if(model.node[colNodeNumber].getMap()==0) 
							v=model.node[colNodeNumber].u;									
						else {
							v=model.node[model.node[colNodeNumber].getMap()].u;
						}

						for(int p=0;p<model.dim;p++)
							if(!model.node[rowNodeNumber].is_U_known(p)){
								rowb=unknown_U_ind[row][p]-1;

								for(int q=0;q<model.dim;q++)
									if(model.node[colNodeNumber].is_U_known(q))
										model.bU.el[rowb]-=Ke[j][k].el[p][q]*v.el[q];

							}


					}

					if(model.node[colNodeNumber].is_U_known()) 
						continue;

					column=model.U_unknownIndex[colNodeNumber]-1;


					if(column>row) continue;

					m=util.search(BKs.row[row].index,nz[row]-1,column);

					if(m<0)
					{	

						BKs.row[row].index[nz[row]]=column;															
						BKs.row[row].el[nz[row]]=Ke[j][k];
						if(massNeeded){
							MKs.row[row].index[nz[row]]=column;	
							MKs.row[row].el[nz[row]]=Me[j][k];
						}

						nz[row]++;

						if(nz[row]==BKs.row[row].nzLength-1){
							BKs.row[row].extend(ext);
							if(massNeeded)
								MKs.row[row].extend(ext);
						}

					}

					else{

						BKs.row[row].addToNz(Ke[j][k],m);
						if(massNeeded)
							MKs.row[row].addToNz(Me[j][k],m);


					}
				}

			}
		}


		BKs.sortAndTrim(nz);


		model.Ks=BKs.spMatForm(model,unknown_U_ind);
		if(massNeeded){

			MKs.sortAndTrim(nz);	
			model.Ms=MKs.spMatForm(model,unknown_U_ind);

		}



	}





	
	public void setNodalMass(Model model){

		System.out.println(" calculating nodal masss... ");

		System.out.println(" Calculating stiffness matrix ...");
		
		double[] nm=new double[model.nElVert];
		
	

		Vect lumpedMass=new Vect(model.numberOfNodes+1);

		for(int i=1;i<=model.numberOfElements;i++){

			if(model.element[i].getRo()==0) {continue;}

			int[] vertNumb=model.element[i].getVertNumb();



				nm=this.calc.nodalMass(model,i);

			for(int j=0;j<model.nElVert;j++){
				
				lumpedMass.el[vertNumb[j]]+=nm[j];
			}
		}

	
		
		for(int i=1;i<=model.numberOfNodes;i++)
			model.node[i].setNodalMass(lumpedMass.el[i]);





	}

	

	public Vect getDeformation( Model model,SpMatSolver solver,int mode){

		solver.terminate(false);
		System.out.println(" Calculating deformation....");


	/*	model.Ks.matForm(true).show();
		
		Eigen eg=new Eigen(model.Ks.matForm(true));
		eg.lam.show();
*/
		Vect bU1=model.bU.add(model.getbUt(mode));

		
		if(model.Ci==null)
			model.Ci=model.Ks.scale(bU1);
		else
			bU1.timesVoid(model.Ci);
		
		Vect u;
		if(model.Ls==null)
		model.Ls=model.Ks.ichol();


		if(bU1.abs().max()!=0){
		if(model.xp==null){
			u=solver.ICCG(model.Ks,model.Ls, bU1,model.errCGmax,model.iterMax);
		}
		else{
			//	u=solver.ICCG(model.Ks,model.Ls, bU1,2e-3,model.iterMax,model.xp);
			u=model.solver.err0ICCG(model.Ks,model.Ls, bU1,model.errCGmax*1e-3,model.iterMax,model.xp);	

		}
		}
		else
			u=new Vect(model.Ks.nRow);

		model.xp=u.deepCopy();

		u.timesVoid(model.Ci);
		

		return u;
	}

	public Vect getVibration( Model model,SpMatSolver solver,int mode,int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using backward Euler two step method....");

		Vect bU1=model.bU.add(model.getbUt(mode));

		double b1=1./Math.pow(model.dt,2);

		double b2=1.0/model.dt;
		
		
		if(iter<2 && !model.loadPrevMech ){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u;
			if(bU1.abs().max()==0)
				u=new Vect( bU1.length); 
			else 
				u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			model.up=u.deepCopy();
			model.upp=u.deepCopy();
			return u;
		}


		if(model.loadPrevMech){


	}
			Vect bp=model.Ms.smul(model.up.times(2).sub(model.upp).times(b1)).add(model.Cs.smul(model.up).times(b2));
			
			bU1=bU1.add(bp);
			
			if(iter==2){
				
				model.Ks=model.Ks.addNew(model.Ms.timesNew(b1)).addNew(model.Cs.timesNew(b2));

			}


		SpMat Ks=model.Ks.deepCopy();

		model.Ci=Ks.scale(bU1);

		Vect u;

		SpMat L=Ks.ichol();

		if(model.xp==null){
			u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);
		}
		else{
			u=model.solver.err0ICCG(Ks,L,bU1,1e-6,model.iterMax,model.xp);	

		}
		

		model.xp=u.deepCopy();

		u.timesVoid(model.Ci);

		if(model.up!=null)
			model.upp=model.up.deepCopy();
		model.up=u.deepCopy();


		return u;
	}



	public Vect getVibrationNewmark( Model model,SpMatSolver solver,int mode,int iter){

		

		solver.terminate(false);
		System.out.println(" Calculating vibration using the Newmark method....");

		double dt=model.dt;

		Vect bU1=model.bU.add(model.getbUt(mode));


		if(iter<2 && !model.loadPrevMech){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);


			SpMat L=Ks.ichol();
			
			Vect u;
			if(bU1.abs().max()==0)
				u=new Vect( bU1.length); 			
			else{
				u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			}
			
			u.timesVoid(model.Ci);	
			
			if(iter==1)
				model.ud=u.sub(model.up).times(1.0/dt);	
			
			model.up=u.deepCopy();
			return u;
		}



		double beta=.25;
		double gama=.5;
		double b1=1./beta/Math.pow(dt,2);

		double b2=-1./beta/dt;

		double b3=1-.5/beta;
		double b4=gama*dt*b1;
		double b5=1+gama*dt*b2;
		double b6=dt*(1+gama*b3-gama);

		SpMat Ks=model.Ks.addNew(model.Ms.timesNew(b1)).addNew(model.Cs.timesNew(b4));

			Vect bp=model.Ms.smul(model.up.times(b1).add(model.ud.times(-b2)).add(model.udd.times(-b3)))
			.add(model.Cs.smul(model.up.times(b4).add(model.ud.times(-b5)).add(model.udd.times(-b6))));

			bU1=bU1.add(bp);


		model.Ci=Ks.scale(bU1);

		Vect u;


		SpMat L=Ks.ichol();

		if(model.xp==null){
			u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);
		}
		else{
			u=model.solver.err0ICCG(Ks,L,bU1,1e-6,model.iterMax,model.xp);	

		}



		model.xp=u.deepCopy();

		u.timesVoid(model.Ci);


		Vect ud1=model.ud.deepCopy();
		Vect udd1=model.udd.deepCopy();
		Vect up1;

		if(model.up==null) up1=new Vect(bU1.length);
		else up1=model.up.deepCopy();
		if(iter>1){
			model.ud=u.sub(up1).times(b4).add(ud1.times(b5)).add(udd1.times(b6));	
			model.udd=u.sub(up1).times(b1).add(ud1.times(b2)).add(udd1.times(b3));
		}

		model.up=u.deepCopy();
	


		return u;

	}
	
	public Vect getVibrationSS22( Model model,SpMatSolver solver,int mode,int iter){


		solver.terminate(false);
		System.out.println(" Calculating vibration using the General SS22 method....");

		double dt=model.dt;

		Vect bU1=model.bU.add(model.getbUt(mode));

		
		if(iter<2){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();
			
			Vect u;
			if(bU1.abs().max()==0)
				u=new Vect( bU1.length); 			
			else
				u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			
			if(iter==1)
			model.ud=u.sub(model.up).times(1.0/dt);	

			model.up=u.deepCopy();
	

			return u;
		}


		double teta1=.5;
		double teta2=.5;
		double b1=dt*teta1;

		double b2=0.5*dt*dt*teta2;

		SpMat Ks=model.Ms.addNew(model.Cs.timesNew(b1)).addNew(model.Ks.timesNew(b2));

			Vect bp=model.Cs.smul(model.ud).add(model.Ks.smul(model.up.add(model.ud.times(b1))));

			bU1=bU1.sub(bp);


		model.Ci=Ks.scale(bU1);

		Vect a2dot;


		SpMat L=Ks.ichol();

		if(model.xp==null){
			a2dot=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);
		}
		else{
			a2dot=model.solver.err0ICCG(Ks,L,bU1,1e-6,model.iterMax,model.xp);	

		}



		model.xp=a2dot.deepCopy();

		a2dot.timesVoid(model.Ci);


		Vect u1=model.up.deepCopy();
		Vect ud1=model.ud.deepCopy();
	
		
		Vect u=u1.add(ud1.times(dt)).add(a2dot.times(0.5*dt*dt));	
		model.ud=ud1.add(a2dot.times(dt));
		
		model.up=u.deepCopy();

		return u;

	}
	
	public Vect getVibrationVerlet( Model model,SpMatSolver solver,int mode,int iter){


		solver.terminate(false);
		System.out.println(" Calculating vibration using the General SS22 method....");

		double dt=model.dt;

		Vect bU1=model.bU.add(model.getbUt(mode));

		
		if(iter<2){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();
			
			Vect u;
			if(bU1.abs().max()==0)
				u=new Vect( bU1.length); 			
			else
				u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			
			if(iter==1)
			model.ud=u.sub(model.up).times(1.0/dt);	

			model.up=u.deepCopy();
	

			return u;
		}




		SpMat Ks=model.Ms.deepCopy();

			Vect bp=model.Ms.addNew(model.Ks.timesNew(-dt)).smul(model.ud).add(model.Ks.smul(model.up.times(-dt)));

			bU1=bU1.sub(bp);


		model.Ci=Ks.scale(bU1);

		Vect v;


		SpMat L=Ks.ichol();

		if(model.xp==null){
			v=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);
		}
		else{
			v=model.solver.err0ICCG(Ks,L,bU1,1e-6,model.iterMax,model.xp);	

		}



		model.xp=v.deepCopy();

		v.timesVoid(model.Ci);


	
		
		Vect u=model.up.add(v.times(dt));	
		
		model.ud=v.deepCopy();
		model.up=u.deepCopy();

		return u;

	}

	public Vect getVibrationCentral( Model model,SpMatSolver solver,int mode, int iter){


		solver.terminate(false);
		System.out.println(" Calculating vibration using the central difference method....");

		double dt=model.dt;


		Vect bU1=model.bU.add(model.getbUt(mode));

		double b1=1./Math.pow(dt,2);

		double b2=.5/dt;


		SpMat Ks;

		if(iter>1){
			
			Ks=model.Ms.timesNew(b1).addNew(model.Cs.timesNew(b2));

			Vect bp=model.Ks.smul(model.up.times(-1)).add(model.Ms.smul(model.up.times(2).sub(model.upp)).times(b1)).add(model.Cs.smul(model.up.times(b2)));

			bU1=bU1.add(bp);

		}
		else{
			Ks=model.Ks.deepCopy();
		}

		model.Ci=Ks.scale(bU1);

		Vect u;


		SpMat L=Ks.ichol();

		if( model.xp==null){
			u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);
		}
		else{
			u=model.solver.err0ICCG(Ks,L,bU1,1e-6,model.iterMax,model.xp);	

		}

		model.xp=u.deepCopy();

		u.timesVoid(model.Ci);

		if(model.up!=null)
			model.upp=model.up.deepCopy();
		model.up=u.deepCopy();

		return u;

	}


	public Vect XEuler( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using the explicit Euler method....");

		double dt=model.dt;




		Vect bU1=model.bU.add(model.getbUt(mode));


		if(iter<1){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);


			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);


			u.timesVoid(model.Ci);	
			model.up=u.deepCopy();
			model.ud=new Vect(u.length);
			return u;
		}

		double b1=1.0/dt;

		Vect ff=new Vect(bU1.length);

		for	(int i=0;i<bU1.length;i++){
			ff.el[i]=bU1.el[i];

		}


		SpMat Ks;


		Ks=model.Ms.timesNew(b1);


		Vect rr=Ks.smul(model.ud).sub(model.Ks.smul(model.up)).sub(model.Cs.smul(model.ud));

		ff=ff.add(rr);


		model.Ci=Ks.scale(ff);


		SpMat L=Ks.ichol();
		Vect v;
		if(iter<2)
		{
			v=solver.ICCG(Ks,L, ff,1e-6,model.iterMax);
		}
		else{
			v=model.solver.err0ICCG(Ks,L,ff,1e-6,model.iterMax,model.xp);	

		}



		model.xp=v.deepCopy();

		v.timesVoid(model.Ci);


		Vect u=model.up.add(model.ud.times(dt));


		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}
	
	public Vect IEuler( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using the implicit Euler method....");

		double dt=model.dt;

		Vect bU1=model.bU.add(model.getbUt(mode));

		if(iter<2){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			if(iter==1)
			model.ud=u.sub(model.up).times(1/dt);
			
			model.up=u.deepCopy();

			return u;
		}


		Vect ff=bU1.times(dt).aug(new Vect(bU1.length));
	
		SpMat Ks=new SpMat(2*model.Ms.nRow,2*model.Ms.nRow);


		SpMat K00=model.Ms.addNew(model.Cs.timesNew(dt));	
				
		SpMat K10=model.Ks.timesNew(dt);
			
		for(int i=0;i<K00.nRow;i++)
			Ks.row[i]=K00.row[i].augh(K10.row[i].deepCopy());
		
		for(int i=0;i<K00.nRow;i++)
			Ks.row[i+K00.nRow]=K10.row[i].augh(model.Ks.row[i].times(-1));
	

		Vect rr1=model.Ms.smul(model.ud);
		Vect rr2=model.Ks.smul(model.up.times(-1));

		Vect rr=rr1.aug(rr2);
		
		
		ff=ff.add(rr);

		
		Vect vu;
		boolean lu=true;
		if(lu){
			if(model.bigMat==null){
				model.bigMat=Ks.matForm(true);

				model.bigMat.lu();
			}
			Mat A=model.bigMat.deepCopy();
			
			vu=new MatSolver().solvelu(A,ff);
		}
			else{
	
			model.Ci=Ks.scale(ff);
	
		//	Ks.diagSym().show();
	
			SpMat L=Ks.ichol();
		
			if(model.xp==null)
			{
				vu=solver.ICCG(Ks,L, ff,1e-6,model.iterMax);
			}
			else{
				vu=model.solver.err0ICCG(Ks,L,ff,1e-6,model.iterMax,model.xp);	
	
			}
	
	
	
			model.xp=vu.deepCopy();
	
			vu.timesVoid(model.Ci);
			}

		Vect v=new Vect(bU1.length);
		Vect u=new Vect(bU1.length);
		
		for(int i=0;i<u.length;i++){
			v.el[i]=vu.el[i];
			u.el[i]=vu.el[i+bU1.length];
		}

		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}
	
	public Vect getVibrationDG( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using the implicit Euler method....");

		double dt=model.dt;

		Vect bU1=model.bU.add(model.getbUt(mode));

		if(iter<2){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			if(iter==1)
			model.ud=u.sub(model.up).times(1/dt);
			
			model.up=u.deepCopy();

			return u;
		}


		Vect ff=bU1.times(dt).aug(new Vect(bU1.length));
	
		SpMat Ks=new SpMat(4*model.Ms.nRow,4*model.Ms.nRow);


		SpMat K00=model.Ms.addNew(model.Cs.timesNew(dt));	
				
		SpMat K10=model.Ks.timesNew(dt);
			
		for(int i=0;i<K00.nRow;i++)
			Ks.row[i]=K00.row[i].augh(K10.row[i].deepCopy());
		
		for(int i=0;i<K00.nRow;i++)
			Ks.row[i+K00.nRow]=K10.row[i].augh(model.Ks.row[i].times(-1));
	

		K00=null;
		K10=null;
		
		SpMat Kdg=new SpMat(2*Ks.nRow,2*Ks.nRow);
		
		Vect rr1=model.Ms.smul(model.ud);
		Vect rr2=model.Ks.smul(model.up.times(-1));

		Vect rr=rr1.aug(rr2);
		
		
		ff=ff.add(rr);

		
		Vect vu;
		boolean lu=true;
		if(lu){
			if(model.bigMat==null){
				model.bigMat=Ks.matForm(true);

				model.bigMat.lu();
			}
			Mat A=model.bigMat.deepCopy();
			
			vu=new MatSolver().solvelu(A,ff);
		}
			else{
	
			model.Ci=Ks.scale(ff);
	
		//	Ks.diagSym().show();
	
			SpMat L=Ks.ichol();
		
			if(model.xp==null)
			{
				vu=solver.ICCG(Ks,L, ff,1e-6,model.iterMax);
			}
			else{
				vu=model.solver.err0ICCG(Ks,L,ff,1e-6,model.iterMax,model.xp);	
	
			}
	
	
	
			model.xp=vu.deepCopy();
	
			vu.timesVoid(model.Ci);
			}

		Vect v=new Vect(bU1.length);
		Vect u=new Vect(bU1.length);
		
		for(int i=0;i<u.length;i++){
			v.el[i]=vu.el[i];
			u.el[i]=vu.el[i+bU1.length];
		}

		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}

	public Vect Runge_Kutta( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using explicit Runga-Kutta method....");

		double dt=model.dt;


		
		Vect bU1=model.bU.add(model.getbUt(mode));

		if(iter<1){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			model.up=u.deepCopy();
			model.ud=new Vect(u.length);
			return u;
		}

		
		Vect k1=ff(model,bU1,model.up,model.ud).times(dt);

		Vect k2=ff(model,bU1,model.up,model.ud.add(k1.times(.5))).times(dt);

		Vect k3=ff(model,bU1,model.up,model.ud.add(k2.times(.5))).times(dt);

		Vect k4=ff(model,bU1,model.up,model.ud.add(k3)).times(dt);


		SpMat Ks;

		Ks=model.Ms.deepCopy();

		Vect b=k1.add(k2.times(2)).add(k3.times(2)).add(k4).times(1.0/6).add(Ks.smul(model.ud));


		model.Ci=Ks.scale(b);


		SpMat L=Ks.ichol();

		Vect v;
		if(model.xp==null)
		{
			v=solver.ICCG(Ks,L, b,1e-6,model.iterMax);
		}
		else{

			v=model.solver.err0ICCG(Ks,L,b,1e-6,model.iterMax,model.xp);		
		}

		model.xp=v.deepCopy();

		v.timesVoid(model.Ci);


		Vect u=model.up.add(v.times(dt));


		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}
	
	public Vect IRunge_Kutta( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using implicit Runga-Kutta method....");

		double dt=model.dt;



		Vect bU1=model.bU.add(model.getbUt(mode));

		if(iter<1){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			model.up=u.deepCopy();
			model.ud=new Vect(u.length);
			return u;
		}

		
		Vect k1=ff(model,bU1,model.up,model.ud).times(dt);

		Vect k2=ff(model,bU1,model.up,model.ud.add(k1.times(.5))).times(dt);

		Vect k3=ff(model,bU1,model.up,model.ud.add(k2.times(.5))).times(dt);

		Vect k4=ff(model,bU1,model.up,model.ud.add(k3)).times(dt);


		SpMat Ks;

		Ks=model.Ms.addNew(model.Ks.timesNew(dt*dt/6)).addNew(model.Cs.timesNew(dt/6));

		Vect b=k2.times(2).add(k3.times(2)).add(k4).times(1.0/6).add(model.Ms.smul(model.ud)).sub(model.Ks.smul(model.up.times(dt/6)));


		model.Ci=Ks.scale(b);


		SpMat L=Ks.ichol();

		Vect v;
		if(model.xp==null)
		{
			v=solver.ICCG(Ks,L, b,1e-6,model.iterMax);
		}
		else{

			v=model.solver.err0ICCG(Ks,L,b,1e-6,model.iterMax,model.xp);		
		}

		model.xp=v.deepCopy();

		v.timesVoid(model.Ci);


		Vect u=model.up.add(v.times(dt));


		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}
	
	
	public Vect leapFrog( Model model,SpMatSolver solver,int mode, int iter){

		solver.terminate(false);
		System.out.println(" Calculating vibration using the leap-frog method....");

		double dt=model.dt;




		Vect bU1=model.bU.add(model.getbUt(mode));


		if(iter<2){

			SpMat  Ks=model.Ks.deepCopy();

			model.Ci=Ks.scale(bU1);

			SpMat L=Ks.ichol();

			Vect u=solver.ICCG(Ks,L, bU1,1e-6,model.iterMax);

			u.timesVoid(model.Ci);	
			model.up=u.deepCopy();
			model.ud=new Vect(u.length);
			return u;
		}


		Vect u=model.up.add(model.ud.times(dt));
		
		Vect ff=new Vect(bU1.length);

		for	(int i=0;i<bU1.length;i++){
			ff.el[i]=bU1.el[i];

		}


		SpMat Ks;


		Ks=model.Ms.deepCopy();


		Vect rr=Ks.smul(model.ud).add(model.Ks.smul(u.times(-dt)));

		ff=ff.add(rr);


		model.Ci=Ks.scale(ff);


		SpMat L=Ks.ichol();
		Vect v;
		if(model.xp==null)
		{
			v=solver.ICCG(Ks,L, ff,1e-6,model.iterMax);
		}
		else{
			v=model.solver.err0ICCG(Ks,L,ff,1e-6,model.iterMax,model.xp);	

		}



		model.xp=v.deepCopy();

		v.timesVoid(model.Ci);


	


		model.up=u.deepCopy();
		model.ud=v.deepCopy();

		return u;

	}

	private Vect ff(Model model,Vect f, Vect u, Vect v){

		Vect v1=f.sub(model.Ks.smul(u).add(model.Cs.smul(v)));
		return v1;
	}


}


