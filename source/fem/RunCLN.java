package fem;

import static java.lang.Math.PI;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Formatter;

import femSolver.CLNStaticMagSolver;
import femSolver.StaticElectricSolver;
import femSolver.StaticLinearMagSolver;
import main.Main;
import math.Complex;
import math.SpMat;
import math.Vect;
import math.VectComp;
import math.util;


public class RunCLN {
	
	String outputFolder;
	int networkType=0;
	int numStages;
	public class ElecMode{
		Vect phiMode;
		Vect aMode;
		
		public ElecMode(int dim){
			phiMode=new Vect(dim);
			aMode=new Vect(dim);
		}
		
		public void Add(ElecMode em){
			phiMode=phiMode.add(em.phiMode);
			aMode=aMode.add(em.aMode);
		}
		
		public void TimesVoid(double a){
			phiMode.timesVoid(a);
			aMode.timesVoid(a);
		}
	}

	public static void main(String[] args){
		new Main();
	}

public void run(Model model, Main main){
		
	numStages=model.nCLNstages;
	
	if(model.nCLNstages<0){
		
		model.nCLNstages=-model.nCLNstages;
		run2(model,main);
		return;
	}else if(model.nCLNstages>100){
		model.nCLNstages-=100;
		networkType=1;
		runFoster(model,main);
		return;	
	}
	
	
	outputFolder=model.resultFolder;
	
		CLNStaticMagSolver magsolver= new CLNStaticMagSolver(model);
		
		model.hasJ =true;
		model.setMagBC();
		
		magsolver.setMagMat(model);
			
		StaticElectricSolver phiSolver= new StaticElectricSolver();
		phiSolver.open_vps=true;
		

		phiSolver.setBoundaryCondition(model);
		phiSolver.setMatrix(model);

		
		int nStages=model.nCLNstages;
		
		double[] losses=new double[1];
		
		Vect resistances=new Vect(nStages);
		Vect inductances=new Vect(nStages);
		
		Vect[] magModes=new Vect[nStages];
		Vect[] elecAModes=new Vect[nStages];
		Vect[] elecPhiModes=new Vect[nStages];

		Vect elecAtemp=new Vect(magsolver.magMat.getnCol());
					
		for(int kstage=0;kstage<nStages;kstage++){
			
	//====== Elec

		if(kstage==0){
			phiSolver.setRHS(model);
		}
		else{
			elecAtemp=magModes[kstage-1].times(-1./inductances.el[kstage-1]);

			
			if (kstage > 1 ||(kstage>0 &&!phiSolver.open_vps)){
				elecAtemp=elecAtemp.add(elecAModes[kstage - 1]);
			}
			
			for(int n=1;n<=model.numberOfNodes;n++) {
				if(model.node[n].isPhiVar()&& model.node[n].isPhiKnown())
				model.node[n].setPhi(0);
			}
						
			model.setSolution(elecAtemp);

			phiSolver.setRHS(model,1.);


		}
		
		
		Vect x=phiSolver.solve(model);
		phiSolver.setSolution(model,x);
		
		 
		if(kstage==0 && phiSolver.open_vps) phiSolver.openVPS(model);
	
		

		if (kstage > 1  ||(kstage>0 &&!phiSolver.open_vps)){
			x=x.add(elecPhiModes[kstage - 1]);
		}
	
		
		util.pr("current: "+x.el[x.length-1]);
		
		elecPhiModes[kstage]=x.deepCopy();
		elecAModes[kstage]=elecAtemp.deepCopy();


		model.setSolution(elecAtemp);

		model.setJStatic();
		model.writeJe(model.resultFolder+"\\Je"+kstage+".txt");
		//model.writePhi(model.resultFolder+"\\phi"+kstage+".txt");
		

		magsolver.setRHS(model,losses);
		
		resistances.el[kstage]=1./losses[0];
		//util.pr("R"+kstage+"="+		resistance.el[kstage]);
		
		//======

		//====== Mag
		magsolver.RHS=magsolver.RHS.times(resistances.el[kstage]);
		
		Vect rhs1=magsolver.RHS.deepCopy();
		
		x=magsolver.solve(model);
		
		util.pr("norma of magSolve deltaA--------------> "+x.norm());
		
	//	if (kstage > 0)
	//		x=x.add(magModes[kstage - 1]);

		
		magModes[kstage]=x.deepCopy();
		
		inductances.el[kstage]=rhs1.dot(x);
		
		//util.pr("L"+kstage+"="+inductances.el[kstage]);
			
		//========
		
		//model.setSolution(x);
		//model.setB();
		//model.writeB(model.resultFolder+"\\B"+kstage+".txt");

		}
		
		for(int kstage=0;kstage<nStages;kstage++){
		util.pr(String.format("%12.5E",resistances.el[kstage]));
		util.pr(String.format("%12.5E",inductances.el[kstage]));

		}
		
		WriteCLN(resistances,inductances);
		//resistance.show();
		//inductances.show();
		//model.writePhi(model.resultFolder+"\\phi.txt");
		
		
		WriteImpedance(resistances,inductances,1e-1,1e3,21);
		
		double f0=1e0;
	//	Complex imp=ObtainImpedance(resistance,inductances,1e2);
		
		Complex vs=new Complex(1,0);
		
		VectComp ee=new VectComp(nStages);
		VectComp hh=new VectComp(nStages);
		
		Complex impedance =SolveCLN(vs,resistances, inductances,ee,hh,f0);
	//	SolveCLNCurrentGiven(vs,resistances, inductances,ee,hh,f0);

		VectComp A=new VectComp(magModes[0].length);
		
		for(int i=0; i<nStages;i++){
			VectComp temp=new VectComp(magModes[i]);
			A=A.add(temp.times(hh.el[i]));
		}
		

		Vect Ar=new Vect(A.length);
		Vect Am=new Vect(A.length);
		
		for(int i=0; i<A.length;i++){
			Ar.el[i]=A.el[i].re;
			Am.el[i]=A.el[i].im;
		}
		
			Ar.times(1e6).hshow();
			Am.times(1e6).hshow();

			model.setSolution(Ar);
			model.setB();
			model.writeB(model.resultFolder+"\\Br.txt");
		
			model.setSolution(magModes[1]);
			model.setB();
			model.writeB(model.resultFolder+"\\Bm.txt");
		
	}




	
	public void run2(Model model, Main main){
		
		
		double w0=2*PI*1e0;
		
		outputFolder=model.resultFolder;
		
			CLNStaticMagSolver magsolver= new CLNStaticMagSolver(model);
			
			model.hasJ =true;
			model.setMagBC();
			
			magsolver.setMagMat(model);
				
			StaticElectricSolver phiSolver= new StaticElectricSolver();
			phiSolver.open_vps=true;
			

			phiSolver.setBoundaryCondition(model);
			phiSolver.setMatrix(model);

			
			int nStages=model.nCLNstages;
			
			double[] losses=new double[1];
			
			//Vect resistance=new Vect(nStages);
			//Vect inductances=new Vect(nStages);
			
			Vect[] magModes=new Vect[nStages];
			Vect[] rhs=new Vect[nStages];
		
			//Vect[] elecAModes=new Vect[nStages];
			Vect[] elecPhiModes=new Vect[nStages];

			//Vect elecAtemp=new Vect(magsolver.magMat.getnCol());
				
			phiSolver.setRHS(model);
	

			
			Vect x=phiSolver.solve(model);
			phiSolver.setSolution(model,x);

				 
			if(phiSolver.open_vps) phiSolver.openVPS(model);
			
			elecPhiModes[0]=x.deepCopy();
			
			magsolver.setRHS(model,losses);
			
			rhs[0]=magsolver.RHS.deepCopy();
			//double rdc=1./losses[0];
			
			//magsolver.RHS.timesVoid(1.*rdc);

		//	util.pr(rdc);
			
		//	resistance.el[0]=1;//1./losses[0];
			
		//	inductances.el[0]=-1;
					
			x=magsolver.solve(model);
			
			VectComp magSol=new VectComp(x.length);

			magSol=magSol.add(new VectComp(x));
			
			magModes[0]=x.deepCopy();
			
			util.pr("norma of magSolve deltaA--------------> "+x.norm());

			
			for(int kstage=1;kstage<nStages;kstage++){
				
	
			//	elecAtemp=magModes[kstage-1].times(-1./inductances.el[kstage-1]);

				
			//	if (kstage > 1 ||(kstage>0 &&!phiSolver.open_vps)){
			//		elecAtemp=elecAtemp.add(elecAModes[kstage - 1]);
			//	}
				
				for(int n=1;n<=model.numberOfNodes;n++) {
					if(model.node[n].isPhiVar()&& model.node[n].isPhiKnown())
					model.node[n].setPhi(0);
				}
							
				model.setSolution(magModes[kstage-1].times(-w0));
				

				phiSolver.setRHS(model,1.);
				
			
			x=phiSolver.solve(model); 		
			
			phiSolver.setSolution(model,x);

		//	if (kstage > 1  ||(kstage>0 &&!phiSolver.open_vps)){
		//		x=x.add(elecPhiModes[kstage - 1]);
		//	}
			
			
			util.pr("current: "+x.el[x.length-1]);
			
			//elecPhiModes[kstage]=x.deepCopy();
			//elecAModes[kstage]=elecAtemp.deepCopy();


			//model.setSolution(x);

			model.setJStatic();
			model.writeJe(model.resultFolder+"\\Je"+kstage+".txt");
			//model.writePhi(model.resultFolder+"\\phi"+kstage+".txt");
			

			magsolver.setRHS(model,losses);

			//magsolver.RHS=magsolver.RHS.add(rhs[kstage-1]);
			//rhs[kstage]=magsolver.RHS.deepCopy();
			//====== Mag
		//	magsolver.RHS=magsolver.RHS.times(resistance.el[kstage]);
			
			Vect rhs1=magsolver.RHS.deepCopy();
			
			x=magsolver.solve(model);
			for(int i=0;i<magSol.length;i++){
				//magSol.el[i].re+=x.el[i];
				magSol.el[i].im+=x.el[i];
			}
			//magSol=magSol.add(new VectComp(x).times(new Complex(0,1)));
			//magSol=new VectComp(x).times(new Complex(0,1));
			
		//	x=x.add(magModes[kstage - 1]);
			util.pr("norma of magSolve deltaA--------------> "+x.norm());
		
			magModes[kstage]=x.deepCopy();
			
		//	inductances.el[kstage]=-1;//rhs1.dot(x);
			
			//util.pr("L"+kstage+"="+inductances.el[kstage]);
				
			//========
			
			//model.setSolution(x);
			//model.setB();
			//model.writeB(model.resultFolder+"\\B"+kstage+".txt");

			}
			
			for(int kstage=0;kstage<nStages;kstage++){
		//	util.pr(String.format("%12.5E",resistance.el[kstage]));
		//	util.pr(String.format("%12.5E",inductances.el[kstage]));

			}
			
			//Vect Ar=magModes[0].deepCopy();
			
			//Ar.times(1e6).hshow();
			//Vect Am=magModes[1].times(1);
			//Am.times(1e6).hshow();

			Vect Ar=new Vect(magSol.length);
			Vect Am=new Vect(magSol.length);
			for(int i=0;i<Ar.length;i++){
				Ar.el[i]=magSol.el[i].re;
				Am.el[i]=magSol.el[i].im;
			}
			
			Ar.times(1e6).hshow();
			//Vect Am=magModes[1].times(1);
			Am.times(1e6).hshow();
			
			model.setSolution(Ar);
			model.setB();
			model.writeB(model.resultFolder+"\\Br.txt");
		
			model.setSolution(Am);
			model.setB();
			model.writeB(model.resultFolder+"\\Bm.txt");
			
		//	WriteCLN(resistance,inductances);
			
		//	WriteImpedance(resistance,inductances,1e1,1e6,21);

			//resistance.show();
			//inductances.show();
			//model.writePhi(model.resultFolder+"\\phi.txt");
		}

	public void runFoster(Model model, Main main){
		
	
	
		outputFolder=model.resultFolder;
		
			CLNStaticMagSolver magsolver= new CLNStaticMagSolver(model);
			
			model.hasJ =true;
			model.setMagBC();
			
			magsolver.setMagMat(model);
				
			StaticElectricSolver phiSolver= new StaticElectricSolver();
			phiSolver.open_vps=true;
			

			phiSolver.setBoundaryCondition(model);
			phiSolver.setMatrix(model);

			
			int nStages=model.nCLNstages;
			
			double[] losses=new double[1];
			
			Vect resistance=new Vect(nStages);
			Vect inductances=new Vect(nStages);
			
			Vect[] magModes=new Vect[nStages];
			Vect[] elecAModes=new Vect[nStages];
			Vect[] elecPhiModes=new Vect[nStages];

			Vect elecAtemp=new Vect(magsolver.magMat.getnCol());
						
			for(int kstage=0;kstage<nStages;kstage++){
				
		//====== Elec

			if(kstage==0){
				phiSolver.setRHS(model);
			}
			else{
				elecAtemp=magModes[kstage-1].times(-1./inductances.el[kstage-1]);

				
				if (kstage > 1 ||(kstage>0 &&!phiSolver.open_vps)){
					elecAtemp=elecAtemp.add(elecAModes[kstage - 1]);
				}
				
				for(int n=1;n<=model.numberOfNodes;n++) {
					if(model.node[n].isPhiVar()&& model.node[n].isPhiKnown())
					model.node[n].setPhi(0);
				}
							
				model.setSolution(elecAtemp);

				phiSolver.setRHS(model,1.);


			}
			
			
			Vect x=phiSolver.solve(model);
			phiSolver.setSolution(model,x);

			
			 
			if(kstage==0 && phiSolver.open_vps) phiSolver.openVPS(model);
		
			

			if (kstage > 1  ||(kstage>0 &&!phiSolver.open_vps)){
				x=x.add(elecPhiModes[kstage - 1]);
			}
		
			
			util.pr("current: "+x.el[x.length-1]);
			
			elecPhiModes[kstage]=x.deepCopy();
			elecAModes[kstage]=elecAtemp.deepCopy();


			model.setSolution(elecAtemp);

			model.setJStatic();
			model.writeJe(model.resultFolder+"\\Je"+kstage+".txt");
			//model.writePhi(model.resultFolder+"\\phi"+kstage+".txt");
			

			magsolver.setRHS(model,losses);
			
			resistance.el[kstage]=1./losses[0];
			//util.pr("R"+kstage+"="+		resistance.el[kstage]);
			
			//======

			//====== Mag
			magsolver.RHS=magsolver.RHS.times(resistance.el[kstage]);
			
			Vect rhs1=magsolver.RHS.deepCopy();
			
			x=magsolver.solve(model);
			
			util.pr("norma of magSolve deltaA--------------> "+x.norm());
			
			if (kstage > 0)
				x=x.add(magModes[kstage - 1]);
	
			
			magModes[kstage]=x.deepCopy();
			
			inductances.el[kstage]=rhs1.dot(x);
			
			//util.pr("L"+kstage+"="+inductances.el[kstage]);
				
			//========
			
			//model.setSolution(x);
			//model.setB();
			//model.writeB(model.resultFolder+"\\B"+kstage+".txt");

			}
			
			for(int kstage=0;kstage<nStages;kstage++){
			util.pr(String.format("%12.5E",resistance.el[kstage]));
			util.pr(String.format("%12.5E",inductances.el[kstage]));

			}
			
			WriteCLN(resistance,inductances);
			//resistance.show();
			//inductances.show();
			//model.writePhi(model.resultFolder+"\\phi.txt");
			
			
			WriteImpedanceFoster(resistance,inductances,1e-1,1e3,21);
			
			double f0=1e2;
		//	Complex imp=ObtainImpedance(resistance,inductances,1e2);
			
			Complex vs=new Complex(1,0);
			
			VectComp ee=new VectComp(nStages);
			VectComp hh=new VectComp(nStages);
			
			//SolveCLN(vs,resistance, inductances,ee,hh,f0);
			SolveCLNCurrentGiven(vs,resistance, inductances,ee,hh,f0);
			
			
			VectComp A=new VectComp(magModes[0].length);
			
			for(int i=0; i<nStages;i++){
				VectComp temp=new VectComp(magModes[i]);
				A=A.add(temp.times(hh.el[i]));
			}
			

			Vect Ar=new Vect(A.length);
			Vect Am=new Vect(A.length);
			
			for(int i=0; i<A.length;i++){
				Ar.el[i]=A.el[i].re;
				Am.el[i]=A.el[i].im;
			}

				Ar.times(1e6).hshow();
				Am.times(1e6).hshow();

			
			
			
		}

	
	private void WriteCLN(Vect res, Vect inds) {

		String type="Cauer";
		if(networkType==1)type="Foster";
		String clnFilePath=outputFolder+"\\cln_out"+type+".txt";
		
		try{
		PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(clnFilePath)));
		
		StringBuilder sbuf = new StringBuilder();
		Formatter fmt = new Formatter(sbuf);

		int numPorts = 1;//cln_generator->GetNumPorts();
		int numStages=res.length;
		
		pw.println("*NO_MODES * NO_PORTS * R_TERMINATION * FOR_SPICE * \n");

		fmt.format("   %d \t %d \t %d \t %d %n", numStages, numPorts, 0, 0);

		pw.print(sbuf.toString());

		pw.println("*** LTspice .cir ************\n");


		pw.println("SPICE\n");



		int cln_id=1;

			for (int i = 0; i < numStages; ++i) {

				double resistance = res.el[i];

				double inductance = inds.el[i];
				
				sbuf = new StringBuilder();
				
				fmt = new Formatter(sbuf);
					if (i == 0){
						fmt.format("    R%d_%d \t p%d \t n%d_%d \t %10.5e%n", cln_id, 2 * i, cln_id, cln_id, i + 1, resistance);

					}
						else{
							fmt.format("    R%d_%d \t n%d_%d \t n%d_%d \t %10.5e%n", cln_id, 2 * i, cln_id, i, cln_id, i + 1, resistance);
						}
				pw.print(sbuf.toString());
					
					sbuf = new StringBuilder();
					fmt = new Formatter(sbuf);
					fmt.format("    L%d_%d \t n%d_%d \t g%d \t %10.5e \t Rser=0%n", cln_id, 2 * i + 1, cln_id, i + 1, cln_id, inductance);
					pw.print(sbuf.toString());
			}

			pw.println("END\n");
			
			util.pr("cln data was written to "+clnFilePath);


			pw.close();

		} catch(IOException e){System.out.println("writing cln file failed.");}

	}
	
	
	
	private void SolveCLNCurrentGiven(Complex current, Vect resistances, Vect inductances,VectComp ee, VectComp hh,double frequency) {

		double omega = 2.*PI*frequency;

		int dim = resistances.length;

		Complex j_omega=new Complex(0., omega);

		Vect R=new Vect(dim);
		VectComp jwL=new VectComp(dim);

		for (int k = 0; k < dim; k++) {
			R.el[k] = resistances.el[k];
			jwL.el[k] = j_omega.times(inductances.el[k]);
		}


		Complex admitance=new Complex(0., 0);
		Complex impedance=new Complex(0., 0);
		
	
		int last = dim - 1;
		impedance = new Complex(R.el[last],0);
		
		impedance=impedance.add(jwL.el[last]);

		for (int k = dim - 2; k >= 0; k--) {
			admitance=new Complex(0., 0);
			if (jwL.el[k].norm()>0)
				admitance=admitance.add(jwL.el[k].inv());
			if (impedance.norm()>0)
				admitance=admitance.add(impedance.inv());

			impedance =  new Complex(R.el[k],0);
			if (admitance.norm()>0)
				impedance=impedance.add(admitance.inv());

		}
		
		
		ee.timesVoid(0);
		hh.timesVoid(0);

		ee.el[0] =  current.times(R.el[0]);
	
		Complex voltage=current.times(impedance);
		
		Complex ik = current;
	
		Complex vk = voltage;

		for (int k = 1; k < dim; k++) {
			vk =vk.sub(ik.times(R.el[k - 1]));

			if (vk.norm()<1e-7) break;

			if (jwL.el[k - 1].im!=0)
				hh.el[k - 1] = vk.times(jwL.el[k - 1].inv());
			else
				hh.el[k - 1] = ik;



			ik=ik.sub(hh.el[k - 1]);
			ee.el[k] =ik.times(R.el[k]);
		}

	
			hh.el[dim - 1] = ee.el[dim - 1].times(1./ R.el[dim - 1]);

		
	}
	private Complex SolveCLN(Complex voltage, Vect resistances, Vect inductances,VectComp ee, VectComp hh,double frequency) {

		double omega = 2.*PI*frequency;

		int dim = resistances.length;

		Complex j_omega=new Complex(0., omega);

		Vect R=new Vect(dim);
		VectComp jwL=new VectComp(dim);

		for (int k = 0; k < dim; k++) {
			R.el[k] = resistances.el[k];
			jwL.el[k] = j_omega.times(inductances.el[k]);
		}


		Complex admitance=new Complex(0., 0);
		Complex impedance=new Complex(0., 0);
		
	
		int last = dim - 1;
		impedance = new Complex(R.el[last],0);
		
		impedance=impedance.add(jwL.el[last]);

		for (int k = dim - 2; k >= 0; k--) {
			admitance=new Complex(0., 0);
			if (jwL.el[k].norm()>0)
				admitance=admitance.add(jwL.el[k].inv());
			if (impedance.norm()>0)
				admitance=admitance.add(impedance.inv());

			impedance =  new Complex(R.el[k],0);
			if (admitance.norm()>0)
				impedance=impedance.add(admitance.inv());

		}
		
		Complex current=voltage.times(impedance.inv());
		
		ee.timesVoid(0);
		hh.timesVoid(0);

		ee.el[0] =  current.times(R.el[0]);
		Complex ik = current;
		Complex vk = voltage;

		for (int k = 1; k < dim; k++) {
			vk =vk.sub(ik.times(R.el[k - 1]));

		//	if (vk.norm()<1e-7) break;

			if (jwL.el[k - 1].im!=0)
				hh.el[k - 1] = vk.times(jwL.el[k - 1].inv());
			else
				hh.el[k - 1] = ik;



			ik=ik.sub(hh.el[k - 1]);
			ee.el[k] =ik.times(R.el[k]);
		}

	
			hh.el[dim - 1] = ee.el[dim - 1].times(1./ R.el[dim - 1]);

		return impedance;
	}
	

	private Complex ObtainImpedance(Vect resistances, Vect inductances, double frequency) {

		
		double omega = 2.*PI*frequency;

		int dim = resistances.length;

		Complex j_omega=new Complex(0., omega);

		Vect R=new Vect(dim);
		VectComp jwL=new VectComp(dim);

		for (int k = 0; k < dim; k++) {
			R.el[k] = resistances.el[k];
			jwL.el[k] = j_omega.times(inductances.el[k]);
		}


		Complex admitance=new Complex(0., 0);
		Complex impedance=new Complex(0., 0);
		
	
		int last = dim - 1;
		impedance = new Complex(R.el[last],0);
		
		impedance=impedance.add(jwL.el[last]);

		for (int k = dim - 2; k >= 0; k--) {
			admitance=new Complex(0., 0);
			if (jwL.el[k].norm()>0)
				admitance=admitance.add(jwL.el[k].inv());
			if (impedance.norm()>0)
				admitance=admitance.add(impedance.inv());

			impedance =  new Complex(R.el[k],0);
			if (admitance.norm()>0)
				impedance=impedance.add(admitance.inv());

		}
		
		return impedance;
	}
	
		
	private void WriteImpedance(Vect resistances, Vect inductances, double f1, double f2,int npoints) {

		
		
		Vect freqs = new Vect(npoints);

		if (npoints == 1){
			freqs.el[0] = f1;
		}

		double factor =  Math.pow(f2 / f1, 1. / (npoints - 1));

		double freq = f1;

	
		for (int i = 0; i < npoints; ++i) {

			freqs.el[i] = freq;

			freq *= factor;
		}
		
		String type="Cauer";
		if(networkType==1)type="Foster";

		String zFilePath=outputFolder+"\\impedance_freq"+type+".txt";

		try{
		PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(zFilePath)));
		
		StringBuilder sbuf = new StringBuilder();
		Formatter fmt = new Formatter(sbuf);

		pw.println("    Port 1 ");
		pw.println();
		pw.println("      Frequency(/s)      R(Ohm)         L(H)");

		
		Complex impedane;

		for (int i = 0; i < npoints; ++i) {

			freq=freqs.el[i];
			
			impedane=ObtainImpedance(resistances, inductances,freq);
			
			sbuf = new StringBuilder();
			
			fmt = new Formatter(sbuf);
			
			fmt.format("   %12.8E  %12.8E  %12.8E %n",freq, impedane.re, impedane.im/(2*PI*freq));

			pw.print(sbuf.toString());
		}
	
			pw.close();

		} catch(IOException e){System.out.println("writing cln file failed.");}

	}
	
	
private void WriteImpedanceFoster(Vect resistances, Vect inductances, double f1, double f2,int npoints) {

		
		
		Vect freqs = new Vect(npoints);

		if (npoints == 1){
			freqs.el[0] = f1;
		}

		double factor =  Math.pow(f2 / f1, 1. / (npoints - 1));

		double freq = f1;

	
		for (int i = 0; i < npoints; ++i) {

			freqs.el[i] = freq;

			freq *= factor;
		}

		
		String zFilePath=outputFolder+"\\impedance_freq_foster.txt";

		try{
		PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(zFilePath)));
		
		StringBuilder sbuf = new StringBuilder();
		Formatter fmt = new Formatter(sbuf);

		pw.println("    Port 1 ");
		pw.println();
		pw.println("      Frequency(/s)      R(Ohm)         L(H)");

		
		Complex admittance=new Complex(0,0);
		Complex impedance=new Complex(0,0);
		for (int i = 0; i < npoints; ++i) {

			freq=freqs.el[i];
			admittance=new Complex(0,0);
			for(int j=0;j<resistances.length;j++){

			admittance=admittance.add(new Complex(resistances.el[j],2*PI*freq*inductances.el[j]).inv());
			}
			
			impedance=admittance.inv();
			
			sbuf = new StringBuilder();
			
			fmt = new Formatter(sbuf);
			
			fmt.format("   %12.8E  %12.8E  %12.8E %n",freq, impedance.re, impedance.im/(2*PI*freq));

			pw.print(sbuf.toString());
		}
	
			pw.close();

		} catch(IOException e){System.out.println("writing cln file failed.");}

	}

}
