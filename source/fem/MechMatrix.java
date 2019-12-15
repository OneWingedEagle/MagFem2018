package fem;


import math.BandMat;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpBlockMat;
import math.SpMat;
import math.SpMatAsym;
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
	
	

	public Vect runContact( Model model,SpMatSolver solver,int mode){

		
		SpMat Ks=null;
		SpMat Ksb=null;
		
		SpMat Kc=null;
		SpMat Kcf=null;
		
		double pf=1e8;
		double pft=1e8;
		
		
		
		Vect bU1=model.bU.add(model.getbUt(mode));
		
		bU1.timesVoid(1e5);

		
		Vect u=new Vect(model.Ks.nRow);
		
		Vect lamN=new Vect(model.Ks.nRow);
		Vect lamT=new Vect(model.Ks.nRow);
		
		Vect aug_N=new Vect(model.Ks.nRow);
		Vect aug_T=new Vect(model.Ks.nRow);
		
		
		Vect Fc=new Vect(model.Ks.nRow);
		Vect Fcf=new Vect(model.Ks.nRow);
		solver.terminate(false);
		System.out.println(" Calculating deformation....");

		int numSn=13;
		int numMed=8;
		int[] slaveNodes=new int[numSn];
		for(int k=0;k<slaveNodes.length;k++)
			slaveNodes[k]=359+k;

		int[][] masterEdges=new int[numMed][2];
		for(int k=0;k<masterEdges.length;k++){
			masterEdges[k][0]=406+k;
			masterEdges[k][1]=masterEdges[k][0]+1;
			//util.hshow(masterEdges[k]);
		}
		
		boolean thick=false;
		
		if(thick){
			 numSn=21;
			 numMed=10;
			 slaveNodes=new int[numSn];
			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=106+k;

			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=127+k;
				masterEdges[k][1]=masterEdges[k][0]+1;
				//util.hshow(masterEdges[k]);	
			}
		}
		
		boolean ipm=false;
		if(ipm){
			 numSn=70;
			 numMed=36;
			 
			// numSn=3;
			// numMed=3;
			 slaveNodes=new int[numSn];
/*			 int kx1=0;
			 slaveNodes[kx1++]=1298;
			 slaveNodes[kx1++]=2408;
			 slaveNodes[kx1++]=1272;*/
	
			 String file="D:\\JavaWorks\\FEM problems\\structural\\2D-contact\\ipm\\cont.txt";
			 int[][] data=model.loader.loadIntArray(file, numSn, 1,1);
			
			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=data[k][0];
		
			data=model.loader.loadIntArray(file, numMed, 2,numSn+2);
		
			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=data[k][0];
				masterEdges[k][1]=data[k][1];
				//util.hshow(masterEdges[k]);	
			}	
		
		}
		boolean punch=true;
		if(punch){
			 numSn=22;
			 numMed=20;
			 
			// numSn=3;
			// numMed=3;
			 slaveNodes=new int[numSn];
/*			 int kx1=0;
			 slaveNodes[kx1++]=1298;
			 slaveNodes[kx1++]=2408;
			 slaveNodes[kx1++]=1272;*/
	
			 String file="D:\\JavaWorks\\FEM problems\\structural\\2D-contact\\punch\\cont.txt";
			 int[][] data=model.loader.loadIntArray(file, numSn, 1,1);
			
			for(int k=0;k<slaveNodes.length;k++)
				slaveNodes[k]=data[k][0];
		
			data=model.loader.loadIntArray(file, numMed, 2,numSn+2);
		
			masterEdges=new int[numMed][2];
			for(int k=0;k<masterEdges.length;k++){
				masterEdges[k][0]=data[k][0];
				masterEdges[k][1]=data[k][1];
				//util.hshow(masterEdges[k]);	
			}	
		
		}
		
		 pf=0;

		for(int i=0;i<slaveNodes.length;i++){
			int sn=slaveNodes[i];
			int index=model.U_unknownIndex[sn]-1;
			if(index<0) continue;
			int xind=2*index;
			int yind=2*index+1;
			double val1=model.Ks.row[xind].el[model.Ks.row[xind].nzLength-1];
			double val2=model.Ks.row[yind].el[model.Ks.row[yind].nzLength-1];
			
			if(val1>pf) pf=val1;
			if(val2>pf) pf=val2;
		}

		pf/=slaveNodes.length;//
		pft=1.*pf;
		
		util.pr("pf :"+pf);
		util.pr("pft :"+pft);

		int[] normalIndex=new int[numSn];
		
		Vect[] normals=new Vect[numMed];

		int dof=model.Ks.nRow;
		

		int nnSize=model.numberOfNodes+1;
		SpMatAsym node_node=new SpMatAsym(nnSize,nnSize);
	//	node_node.lower=false;

		
		//node_node.shownzA();
		
		//node_node.matForm().plot();
		
		int dim=model.dim;

		SpMatAsym Gc=new SpMatAsym(dof,dof);
		SpMatAsym Gcf=new SpMatAsym(dof,dof);
		SpMatAsym Gct=null;
		SpMatAsym Gcft=null;
		int itmax=10;
		int nr_itmax=10;
		int nLoads=10;
		
		Vect err=new Vect(itmax);
		
		Vect errf=new Vect(itmax);
		
		Vect nr_err=new Vect(nLoads*itmax*nr_itmax);
		
		int totalNRIter=0;
Vect load=bU1.deepCopy();
for(int load_iter=0; load_iter<nLoads; load_iter++){
	
double factor=(load_iter+1.)/nLoads;
	bU1=load.times(factor);

for(int cont_iter=0; cont_iter<itmax; cont_iter++){
	
util.pr("cont_iter: "+cont_iter);

	for(int nr_iter=0; nr_iter<nr_itmax; nr_iter++){	
		util.pr("nr_iter: "+nr_iter);

		int numContacting=0;
		
		for(int i=0;i<slaveNodes.length;i++){
			int sn=slaveNodes[i];
			Vect v=model.node[sn].getCoord().add(model.node[sn].u);
			for(int k=0;k<masterEdges.length;k++){
				int mn1=masterEdges[k][0];
				int mn2=masterEdges[k][1];
				Vect v1=model.node[mn1].getCoord().add(model.node[mn1].u);
				Vect v2=model.node[mn2].getCoord().add(model.node[mn2].u);
				Vect v12=v2.sub(v1);
				//v12.hshow();
				Vect v1v=v.sub(v1);
				Vect v2v=v.sub(v2);
				
				Vect edgeDir=v12.normalized();
		
				double dot1=v1v.dot(v12);
				double dot2=v2v.dot(v12);
				if(dot1*dot2>0) continue;
				
		
				Vect normal=new Vect(-v12.el[1],v12.el[0]).normalized();
				normals[k]=normal;
				//if(v1.el[1]<.0201) 
					normals[k].timesVoid(-1);
			//	normals[k].hshow();
				
				double pen=v1v.dot(normal);
				
				if(pen>0) continue;
				
				double edgeLength=v12.norm();

				if(pen<-10*edgeLength) continue;
				
				normalIndex[i]=k;
	
				double beta=0;
				double alpha=0;
				double v1n=v1v.norm();

				if(v1n==0){
					 beta=0;
					 alpha=1;
				}
				///if(sn>=120) util.pr(sn+"   "+v1n);
				else{
				beta=v1v.dot(edgeDir)/edgeLength;			
				alpha=1-beta;
				}
				//util.pr((sn)+" -edgeLength--  "+ edgeLength);
				//util.pr((sn)+" ---  "+ mn1+"  "+mn2);
				//util.pr((sn)+" ---  "+ alpha+"  "+beta);
				Vect relation=new Vect(nnSize);
				relation.el[mn1]=alpha;
				relation.el[mn2]=beta;
				
				node_node.row[sn]=new SpVect(nnSize,2);
				node_node.row[sn].index[0]=mn1;
				node_node.row[sn].index[1]=mn2;
				node_node.row[sn].el[0]=alpha;
				node_node.row[sn].el[1]=beta;
		
				numContacting++;
			//	break;//
			}
			
		}
		
		util.pr("numContacting:   "+numContacting);

		for(int i=0;i<slaveNodes.length;i++){
			int sn=slaveNodes[i];

			if(node_node.row[sn].nzLength>0){
				
				int index=model.U_unknownIndex[sn]-1;
				if(index<0) continue;
				
				
				int com_index=dim*index;

				
				Vect normal=normals[normalIndex[i]];
				
				int mn1=masterEdges[normalIndex[i]][0];
				int mn2=masterEdges[normalIndex[i]][1];
				
				
				int index1=model.U_unknownIndex[mn1]-1;
				int index2=model.U_unknownIndex[mn2]-1;
				

				int com_index1=dim*index1;
				int com_index2=dim*index2;
				
				double alpha=node_node.row[sn].el[0];
				double beta=node_node.row[sn].el[1];
				

				Gc.row[com_index]=new SpVect(nnSize,6);
				int kx=0;
				Gc.row[com_index].index[kx]=com_index;
				Gc.row[com_index].el[kx++]=normal.el[0];
				Gc.row[com_index].index[kx]=com_index+1;
				Gc.row[com_index].el[kx++]=normal.el[1];
		
				if(com_index1>=0){
				Gc.row[com_index].index[kx]=com_index1;
				Gc.row[com_index].el[kx++]=-alpha*normal.el[0];
				Gc.row[com_index].index[kx]=com_index1+1;
				Gc.row[com_index].el[kx++]=-alpha*normal.el[1];
				}
				if(com_index2>=0){
				Gc.row[com_index].index[kx]=com_index2;
				Gc.row[com_index].el[kx++]=-beta*normal.el[0];
				Gc.row[com_index].index[kx]=com_index2+1;
				Gc.row[com_index].el[kx++]=-beta*normal.el[1];
				}
				Gc.row[com_index].sortAndTrim(kx);;

				
			//	util.pr((sn)+"    "+(com_index)+"  ===== "+ (com_index1+1)+"  "+(com_index2+1));
			//	util.pr((sn)+"    "+(com_index)+"  ===== "+ alpha+"  "+beta);

				
				Vect tang=new Vect(-normal.el[1],normal.el[0]);
				Gcf.row[com_index]=new SpVect(nnSize,6);
				 kx=0;
				Gcf.row[com_index].index[kx]=com_index;
				Gcf.row[com_index].el[kx++]=tang.el[0];
				Gcf.row[com_index].index[kx]=com_index+1;
				Gcf.row[com_index].el[kx++]=tang.el[1];
		
				if(com_index1>=0){
				Gcf.row[com_index].index[kx]=com_index1;
				Gcf.row[com_index].el[kx++]=-alpha*tang.el[0];
				Gcf.row[com_index].index[kx]=com_index1+1;
				Gcf.row[com_index].el[kx++]=-alpha*tang.el[1];
				}
				if(com_index2>=0){
				Gcf.row[com_index].index[kx]=com_index2;
				Gcf.row[com_index].el[kx++]=-beta*tang.el[0];
				Gcf.row[com_index].index[kx]=com_index2+1;
				Gcf.row[com_index].el[kx++]=-beta*tang.el[1];
				}
				
				Gcf.row[com_index].sortAndTrim(kx);;

	
				//if(mn1==408) util.pr((com_index1+1)+" ---  "+ alpha+"  "+beta);
			//	if(mn2==408) util.pr((com_index2+1)+" ---  "+ alpha+"  "+beta);
				//util.pr(sn+"   "+index);
				//Gc.row[index].showr();
				//Gc.row[index].shownz();
			}
		}
	
		if(numContacting!=0){
			
		//Gc.shownzA();
	//	Gct=new SpMatAsym(Gc.matForm().transp());
		Gct=Gc.transpose(500);
		
		//Gct.shownzA();
	
		Kc=new SpMat(dof,dof); // Gct*Gc
		
		for(int i=0;i<Gct.nRow;i++){
			if(Gct.row[i].nzLength>0){
				SpVect spv=new SpVect(dof,100);
			
				int kx=0;
				for(int j=0;j<=i;j++){
					if(Gct.row[j].nzLength>0){
						
						double dot=Gct.row[i].dot(Gct.row[j]);
						if(dot==0) continue;
					//	util.pr(i+" ---- "+j);
						spv.index[kx]=j;
						spv.el[kx++]=dot;
					}
				}
				spv.trim(kx);
				Kc.row[i]=spv.deepCopy();
				
			}
		}
		
		Kc.times(pf);
		
		
		Gcft=Gcf.transpose(500);
		
Kcf=new SpMat(dof,dof); // Gct*Gc
		
		for(int i=0;i<Gcft.nRow;i++){
			if(Gcft.row[i].nzLength>0){
				SpVect spv=new SpVect(dof,100);
			
				int kx=0;
				for(int j=0;j<=i;j++){
					if(Gcft.row[j].nzLength>0){
						
						double dot=Gcft.row[i].dot(Gcft.row[j]);
						if(dot==0) continue;
					//	util.pr(i+" ---- "+j);
						spv.index[kx]=j;
						spv.el[kx++]=dot;
					}
				}
				spv.trim(kx);
				Kcf.row[i]=spv.deepCopy();
				
			}
		}
		
		Kcf.times(pft);
	//	Kc.shownzA();
	//	Kc.plot();
		
		//Mat M=model.Ks.matForm();
		
	//	M=M.add(Kc.matForm());
		
		//Ks= new SpMat(M);

	
		//if(Ksb==null){
		 Ks=model.Ks.addGeneral(Kc.addGeneral(Kcf));
			Ksb=Ks.deepCopy();

		//}
	//	 else
		//	Ks=Ksb.deepCopy();
		 
		
		// Ks.plot();
	

		}
		else
		Ks=model.Ks.deepCopy();
	
		//Ks.plot();
		
	///	util.wait(10000);
		//model.Ks.matForm(true).show();
		
	//	model.Ks.showcl();
	
		boolean solve=true;
		
		
		if(solve){

		
		Vect Fint=Ks.smul(u);


		Vect dF=bU1.sub(Fint);
		
		//Gct.shownzA();
				
		if(Kc!=null){


		//Fc=Kc.smul(u).times(.5);
		
		//aug_N=aug_N.add(Fc);
		
		//Fcf=Kcf.smul(u).times(.5);
		
		//aug_T=aug_T.add(Fcf);


		dF=dF.sub(aug_N).sub(aug_T);



		}

	
		
		double er=dF.norm()/bU1.norm();
		nr_err.el[totalNRIter]=er;
		
		if(er<1e-3) break;
		
		totalNRIter++;


		model.Ci=Ks.scale(dF);
	
		
		Vect du;
		//if(model.Ls==null)
		model.Ls=Ks.ichol();


		if(dF.abs().max()!=0){
		//if(model.xp==null){
			du=solver.ICCG(Ks,model.Ls, dF,model.errCGmax,model.iterMax);
	//	}
		//else{
			//	u=solver.ICCG(model.Ks,model.Ls, bU1,2e-3,model.iterMax,model.xp);
		//	du=model.solver.err0ICCG(Ks,model.Ls, dF,model.errCGmax*1e-3,model.iterMax,model.xp);	

		//}
		}
		else{
			util.pr("Solution is zero!");
			du=new Vect(Ks.nRow);
		}

		model.xp=du.deepCopy();

		du.timesVoid(model.Ci);
		
		u=u.add(du);
		
		model.setU(u);
		
			
		}
		
				
	}
		Vect gap=Gc.mul(u);
		

//		lamN=lamN.add(Fc);
		//lamN=lamN.add(Kc.smul(u));

		err.el[cont_iter]=gap.abs().max();
		
		
		
		//lamT=lamT.add(Kcf.smul(u));
		
		Vect sld=Gcf.mul(u);
		errf.el[cont_iter]=sld.abs().max();
		
		if(err.el[cont_iter]<1e-6 && (errf.el[cont_iter]<1e-4)) break;
		
		if(Kc!=null){
			
			Vect g=Gc.mul(u).times(pf);
			Vect s=Gcf.mul(u).times(pft);
			lamN=lamN.add(g);
			lamT=lamT.add(s);
			aug_N=Gct.mul(lamN);;
			aug_T=Gcft.mul(lamT);;

		}


	
}

	}


util.pr("NR error");

int nn=0;
for(int k=0;k<nr_err.length;k++)
	if(nr_err.el[k]>0) nn++;

Vect nr_distilled=new Vect(nn);

int kx=0;
for(int k=0;k<nr_err.length;k++)
	if(nr_err.el[k]>0) nr_distilled.el[kx++]=nr_err.el[k];

nr_distilled.show();
util.plot(nr_distilled);
//u=aug_N.add(aug_T);


util.pr("Gap[micon] vs aug_iter");
err.times(1e6).show();
//util.plot(err);
	
util.pr("slide[micon] vs aug_iter");
errf.times(1e6).show();
//util.plot(errf);


return u;

	}
	
}


