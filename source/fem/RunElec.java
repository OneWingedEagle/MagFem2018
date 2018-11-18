package fem;

import static java.lang.Math.PI;

import java.io.File;

import java.text.DecimalFormat;

import femSolver.StaticElectricSolver;
import main.Main;
import math.Vect;
import math.util;

public class RunElec {

	private DecimalFormat formatter=new DecimalFormat("0.00");

	public static void main(String[] args){
		new Main();
	}

	public void run(Model model, Main main){
		/*
		StaticElectricSolver esolver=new StaticElectricSolver();
		
		esolver.setBoundaryCondition(model);
		esolver.setMatrix(model);
		esolver.setRHS(model);
		esolver.solve(model, 1);*/
		
		for(int ir=1;ir<=1;ir++){
		PhiCoil coil=new PhiCoil(ir);
		double sig=model.region[ir].getSigma().el[0];
		coil.setSigma(sig);
		
		coil.setBoundaryCondition(model);
		coil.setMatrix(model);
		coil.setRHS(model);

		Vect x=coil.solve(model);
		
		model.setJStatic();
		}
		
		model.writePhi(model.resultFolder+"\\phi.txt");
		model.writeJe(model.resultFolder+"\\Je.txt");
	}

}
