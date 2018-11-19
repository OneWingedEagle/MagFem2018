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
import math.SpMat;
import math.Vect;
import math.util;

public class RunCLN {
	
	String outputFolder;

	public static void main(String[] args){
		new Main();
	}
	
	
public void run(Model model, Main main){
		
	
	outputFolder=model.resultFolder;
	
		CLNStaticMagSolver magsolver= new CLNStaticMagSolver(model);
		
		model.hasJ =true;
		model.setMagBC();
		
		magsolver.setMagMat(model);
			
		StaticElectricSolver phiSolver= new StaticElectricSolver();
		phiSolver.open_vps=true;
		

		phiSolver.setBoundaryCondition(model);
		phiSolver.setMatrix(model);

		
		int nStages=3;
		
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
	}



	public void runbb(Model model, Main main){
		
		CLNStaticMagSolver magsolver= new CLNStaticMagSolver(model);
		
		model.hasJ =true;
		model.setMagBC();
		
		//boolean openVPS=false; // not working with false
		
		magsolver.setMagMat(model);
		
		boolean together=true;
		
		StaticElectricSolver phiSolver= new StaticElectricSolver();
		phiSolver.open_vps=true;
		
		if(together){
		phiSolver.setBoundaryCondition(model);
		phiSolver.setMatrix(model);
		}
		else{
		for(int ic=0;ic<model.phiCoils.length;ic++){
		PhiCoil coil=model.phiCoils[ic];
		double sig=model.region[1].getSigma().el[0];
		coil.setSigma(sig);
		
		
		coil.setBoundaryCondition(model);
		coil.setMatrix(model);
		}
		}
		
		int nStages=3;
		
		double[] losses=new double[1];
		
		Vect resistance=new Vect(nStages);
		Vect inductances=new Vect(nStages);
		
		Vect[] magModes=new Vect[nStages];
		Vect[] elecAModes=new Vect[nStages];
		Vect[] elecPhiModes=new Vect[nStages];

		Vect elecAtemp=new Vect(magsolver.magMat.getnCol());
		
		PhiCoil coil;
		
		coil=model.phiCoils[0];
		
		for(int kstage=0;kstage<nStages;kstage++){
			
	//====== Elec

		if(kstage==0){
		//
			if(together)
			phiSolver.setRHS(model);
			else
				coil.setRHS(model);
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
			if(together)
				phiSolver.setRHS(model,1.);
				else
					coil.setRHS(model,1.);

		}
		
		
		Vect x=null;
		if(together){
			 x=phiSolver.solve(model);
			 if(kstage==0 && phiSolver.open_vps) phiSolver.openVPS(model);
		}
		else x=coil.solve(model);
		

		if (kstage > 1  ||(kstage>0 &&!phiSolver.open_vps)){
			x=x.add(elecPhiModes[kstage - 1]);
		}
		
		
		util.pr("current: "+x.el[x.length-1]);
		
		elecPhiModes[kstage]=x.deepCopy();
		elecAModes[kstage]=elecAtemp.deepCopy();


		model.setSolution(elecAtemp);
	//	x.show();
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
		//resistance.show();
		//inductances.show();
		//model.writePhi(model.resultFolder+"\\phi.txt");
	}

	
	private void WriteCLN(Vect res, Vect inds) {

		String clnFilePath=outputFolder+"\\cln_out.txt";
		
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
	

}
