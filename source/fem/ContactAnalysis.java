

package fem;

import static java.lang.Math.PI;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import io.Loader;
import main.Main;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatAsym;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;


public class ContactAnalysis {
	

	Contact contact;

	private SpMat Ks = null;


	private SpMat Kc = null;
	private SpMat Kcf = null;

	Mat KK = null;
	MatSolver direct_slv = null;

	private Vect aug_N;
	private Vect aug_T;

	private Vect ref_stick;

	private Vect Fc;
	private Vect Fcf;

	private double penalMax;

	int[][] u_index;
	int[] u_index_inv;
	int totalnumContactingNodes;


	Model model;
	double fr = 1e-2;
//	double pft = 1e8;
	double margin = 0.02;

	int n_modifNR0 = 5;
	int aug_itmax0 = 3;
	int nr_itmax0 = 10;
	int nLoads0 = 1;
	double nr_tol0=1e-3;
	double aug_tol0=1e-6;
	
	double nr_tol;
	double aug_tol;
	double mnr_tol;
	
	int n_modifNR;
	int aug_itmax;
	int nr_itmax;
	int nLoads;
	
	Vect disp, rhs;
	boolean applyNodal = true;
	boolean initialized = false;


	boolean aug_normal = true;
	boolean aug_tang = true;

	double extention_fact = .01;
	double clrFact = 1e-10;
	double aug_disp_tol = 1e-4;
	double gap_tol = 1e-4;
	double reduct = .0;


	double adh = 0e3;
	double adhf = 0e8;

	double relax =.5;
	
	Vect top;
	SpMat K_hat;
	Vect rhs_hat;
	int totalNRIter = 0;
	boolean plot_radial=false;
	
	Vect gap_err,slide_err, aug_disp_err,nr_err;
	int nr_it[];


	public Vect solve(Model model, SpMatSolver solver,SpMat Khat1,Vect bhat1,int step) {
		
		this.K_hat=Khat1.deepCopy();
		this.rhs_hat=bhat1.deepCopy();

		if(model.dim==3) 	extention_fact = .02;

		
		if(step==model.nTsteps-1) plot_radial=true;
		
		solver.terminate(false);
		this.model = model;

		direct_slv = new MatSolver();



		fr = .01;



		applyNodal = true;

		
		KK = null;

		boolean direct = true;
	


			if(!initialized){
				
				disp = new Vect(model.Ks.nRow);
				if (disp.length > 10000 || n_modifNR==0)
					direct = false;
				
			initialize();
			initialized=true;
			
			calcPenaltyFactor();
			}

			rhs =rhs_hat.deepCopy();

			System.out.println(" Contact analysis....");

		//	int[] mm;
			int nout = 0;
			int[] xr_nids1=new int[contact.slaveNodes[0].length];

			for (int k = 0; k <contact.slaveNodes[0].length; k++) {
				Node snode=contact.slaveNodes[0][k];
				if(model.dim==2 ||snode.getCoord(2)<1e-6){
					xr_nids1[nout]=snode.id;
					nout++;
				}
			}

			util.pr(nout+"  <==== ");
			Vect xr = new Vect(nout);

	
			if (plot_radial)
				for (int k = 0; k <xr.length; k++) {
					int sn=xr_nids1[k];
					Node snode=model.node[sn];
					if(model.dim==2 ||snode.getCoord(2)<1e-4){
						xr.el[k] = snode.getCoord(0);

					}
				}
			//util.hshow(xr_nids1);
		//xr.hshow();
			int[] indx=xr.bubble();
			int[] xr_nids=new int[nout];
			for (int k = 0; k < xr.length; k++)
				xr_nids[k]=	xr_nids1[indx[k]];
			//xr.hshow();
			//util.hshow(xr_nids);


			Vect dF = null;

			Vect[] urs = new Vect[aug_itmax];
			for (int i = 0; i < aug_itmax; i++)
				urs[i] = new Vect(xr.length);


		
			Vect load = rhs.deepCopy();
			for (int load_iter = 0; load_iter < nLoads; load_iter++) {

				util.pr("load_iter: " + load_iter);

				double factor = (load_iter + 1.) / nLoads;
				rhs = load.times(factor);
				
				for (int aug_iter = 0; aug_iter < aug_itmax; aug_iter++) {
					

					util.pr("aug_iter: " + aug_iter);

					Vect uaug = disp.deepCopy();

					double er = 1;
					double disp_err = 1;

					boolean zigzag = false;

					for (int nr_iter = 0; nr_iter < nr_itmax; nr_iter++) {
						
						assembleContactMatrices();
						
						//util.pr(" here 666");
						addMatrices();
						
						
						boolean mnr_converged =handleModifiedMNR(nr_iter,solver,direct);


						if(mnr_converged) break;
					

						dF = this.calcResidual(Ks, disp, rhs,aug_iter);
		
						er = dF.norm() / rhs.norm();
						nr_err.el[totalNRIter] = er;
						nr_it[totalNRIter] = totalNRIter;
						

						Vect du = solveLinear(solver, Ks.deepCopy(), dF);

						disp = disp.add(du);

						model.setU(disp);

					//	model.writeNodalField( model.eddyFolder+"\\disp_nr"+totalNRIter+".txt",-1);
					
						if (disp.norm() > 0)
							disp_err = du.norm() / disp.norm();

						util.pr("nr_iter: " + nr_iter + "           nr_err: " + er + "       disp_err: " + disp_err);


						totalNRIter++;
						if (er < nr_tol) {
							break;
						}
							
					}

					double dip_err = disp.sub(uaug).norm() / disp.norm();

					aug_disp_err.el[aug_iter] = dip_err;

					if (dip_err < aug_tol)
						break;


					updateLambdas();

					double[] aug_errs=getGapAndSlideErrors();

					checkStickSlip();

					updateAugLoads();



					gap_err.el[aug_iter] = aug_errs[0];
					slide_err.el[aug_iter] = aug_errs[1];



					if(gap_err.el[aug_iter]<aug_tol && (!contact.frictional ||aug_disp_err.el[aug_iter]<aug_tol))
					{
						break;
					}

					
				gap_err.el[aug_iter] = aug_errs[0];
				slide_err.el[aug_iter] = aug_errs[1];

					// xr.show();
					if (plot_radial)
						plotRadialDisp( xr, urs, xr_nids,  aug_iter);
						
						 
					 if(gap_err.el[aug_iter]<aug_tol && (!contact.frictional ||aug_disp_err.el[aug_iter]<aug_tol))
					 {
				//	 break;
					 }

				}
				
			//	model.writeNodalField( model.eddyFolder+"\\disp_nr"+load_iter+".txt",-1);
			}


			
			if(step==model.nTsteps-1){
				printoutErrors();
			
			if (plot_radial) {
				double[][] data = new double[xr.length][aug_itmax + 1];
				for (int k = 0; k < xr.length; k++) {

					data[k][0] = xr.el[k];

					for (int i = 0; i < aug_itmax; i++)
						data[k][i + 1] = urs[i].el[k];

				}
				// util.show(data);
				if(data.length>0)
				util.plotBunch(data);
				// util.plot(x,ur);
			}
			


		}
	
		writeContactForce();


		model.setU(disp);

		return disp;

	}

	private void assembleContactMatrices() {

		int dof = model.Ks.nRow;

		if (model.dim == 2) {
			obtain_node_node();
			assembleConstraintMats();
		} else {
			obtain_node_node3D();
		//	util.pr(" here 444");
			assembleConstraintMats3D();
			//node_node_mat[contId].shownzA();
		///	contact.constraint_matrix_N[contId].shownzA();
		}

		util.pr("totalnumContactingNodes: " + totalnumContactingNodes);

		

		//Kadh = new SpMat(dof, dof);

		Kc = new SpMat(dof, dof); // Gct*Gc
		
		if (totalnumContactingNodes != 0) {
///contact.constraint_matrix_N[contId].shownzA();
	
			for (int contId = 0; contId < contact.numContacts; contId++){
		//	util.pr(" here 1111");
			contact.constraint_matrix_N_trp[contId] = contact.constraint_matrix_N[contId].transpose(100);

			SpMat Kc1=new SpMat(dof, dof);
			SpMat Kadh1=new SpMat(dof, dof);
			
			for (int i = 0; i < contact.constraint_matrix_N_trp[contId].nRow; i++) {
				if (contact.constraint_matrix_N_trp[contId].row[i].nzLength > 0) {
					SpVect spv =null;// new SpVect(dof, 100);
					SpVect spvx =null; //new SpVect(dof, 100);

					boolean spv1_filled=false;
					SpVect spv1 =null;// contact.constraint_matrix_N_trp[contId].row[i].deepCopy();

					SpVect spv2 =null;// contact.constraint_matrix_N_trp[contId].row[i].deepCopy();
				

					int kx = 0;

					for (int j = 0; j <= i; j++) {
						if (contact.constraint_matrix_N_trp[contId].row[j].nzLength > 0) {
							
							if(!spv1_filled){
								 spv = new SpVect(dof, 100);
								 spvx = new SpVect(dof, 100);
								 spv1 = contact.constraint_matrix_N_trp[contId].row[i].deepCopy();

								spv2 =contact.constraint_matrix_N_trp[contId].row[i].deepCopy();
								for (int k = 0; k < spv1.nzLength; k++) {
									int ind = spv1.index[k];
								//	util.pr(ind+"  / "+weights.length);
									spv1.el[k] *= contact.weights[contId].el[ind];
								}

								spv1_filled=true;
							}
							double dot = spv1.dot(contact.constraint_matrix_N_trp[contId].row[j]);
							double dotx = spv2.dot(contact.constraint_matrix_N_trp[contId].row[j]);

							if (dot == 0)
								continue;

							spvx.index[kx] = j;
							spvx.el[kx] = dotx;

							spv.index[kx] = j;
							spv.el[kx++] = dot;

						}

					}
					if(spv!=null){
					spv.trim(kx);
					Kc1.row[i] = spv.deepCopy();
					

					spvx.trim(kx);
					Kadh1.row[i] = spvx.times(adh);
					}

				}
			}
			
			if(Kc1.size()>0){
				Kc=Kc.addGeneral(Kc1);
			}
			
			if(Kadh1.size()>0){
				Kc=Kc.addGeneral(Kadh1);
			}
			
			}
			
		//	util.pr(" here 222");

		//	Kc.times(pf);

			
			Kcf = new SpMat(dof, dof); // Gct*Gc

			// === adhisve tang
		//	Kadhf = new SpMat(dof, dof);

			for (int contId = 0; contId < contact.numContacts; contId++){
				
				
		//	util.pr(" here 1111");
			contact.constraint_matrix_T_trp[contId] = contact.constraint_matrix_T[contId].transpose(100);
			
			contact.constraint_matrix_STK_trp[contId] = contact.constraint_matrix_STK[contId].transpose(100);



			SpMat Kcf1=new SpMat(dof, dof);
	


			for (int i = 0; i < contact.constraint_matrix_STK_trp[contId].nRow; i++) {
				if (contact.constraint_matrix_STK_trp[contId].row[i].nzLength > 0) {
					

					
					SpVect spv = new SpVect(dof, 100);
					SpVect spv1 = contact.constraint_matrix_STK_trp[contId].row[i].deepCopy();

					for (int k = 0; k < spv1.nzLength; k++) {
						int ind0 = spv1.index[k];

						int ind=u_index_inv[ind0];
						spv1.el[k] *= contact.weights[contId].el[ind];
					}

					int kx = 0;
					for (int j = 0; j <= i; j++) {
						if (contact.constraint_matrix_STK_trp[contId].row[j].nzLength > 0) {

							double dot = spv1.dot(contact.constraint_matrix_STK_trp[contId].row[j]);

							if (dot == 0)
								continue;
							spv.index[kx] = j;
							spv.el[kx++] = dot;
						}
					}
					spv.trim(kx);
					Kcf1.row[i] = spv.deepCopy();

				}
			}
			
			if(Kcf1.size()>0){
				Kcf=Kcf.addGeneral(Kcf1);
			}
			

			}

			Kcf.times(fr);

		}
		
	//	util.pr(" KcFFFFFFFFFFFFFFFFFFFFFFFFFF   "+Kcf.norm());
		
		//util.pr(" here 3333");
	}

	private void addMatrices() {

		if (Kc != null) {
			if(Kcf!=null)
				Kc=Kc.addGeneral(Kcf);
		///	Kcf.shownzA();

			Ks = K_hat.addGeneral(Kc);
			
		//	Ks = Ks.addGeneral(Kadh.addGeneral(Kadhf));
		} else {
			Ks = K_hat.deepCopy();
		}
		
	}

	private double modifiedNR(SpMatSolver solver, Vect bU1, boolean direct, double tol) {

		double er1 = 1;
		Vect uc[] =new Vect[n_modifNR];
		Vect errs=new Vect().ones(n_modifNR).times(-1);
		for (int sb = 0; sb < n_modifNR; sb++) {


			 
			if(model.dim==2){
				
				  if(sb%1==0){ 
					 assembleContactMatrices(); 
					 addMatrices();
		 
				  }else
					  updateGap(true);
			}
			else
			{
				obtain_node_node3D();
				assembleConstraintMats3D();
		
			}
	

			Vect dF1 = calcResidual(Ks, disp, bU1,0);

			er1 = dF1.norm() / bU1.norm();

			errs.el[sb]=er1;
			uc[sb]=disp.deepCopy();
			
			util.pr("                                                      nr_iter_sub: " + sb
					+ "           nr_err_sub: " + er1);

			if (sb > 1 && er1 > 100)
				break;

			if (er1 < tol) {
				break;
			}

			Vect du1 = null;

			if (direct) {
				du1 = direct_slv.solvelu(KK, dF1);
			} else
				du1 = solveLinear(solver, Ks.deepCopy(), dF1.deepCopy());


			disp = disp.add(du1);

			model.setU(disp);

		}
		
		int min_index=0;
		double min_err=1e10;
		for(int i=0;i<errs.length;i++){
			if(errs.el[i]>=0 && errs.el[i]<min_err) {
				min_err=errs.el[i];
				min_index=i;
			}
		}

		//if (er1 > tol) {
			disp = uc[min_index].deepCopy();
			model.setU(disp);
		//}
			
			return er1;

	}

	private void obtain_node_node() {



		// weights.zero();
		// weightsf.zero();


		totalnumContactingNodes = 0;

		int nnSize = model.numberOfNodes + 1;
		for (int contId = 0; contId < contact.numContacts; contId++) {
			
			contact.numContactingNodes[contId]=0;
			
			contact.gap[contId].zero();
			contact.slide[contId].zero();
			
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {
				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;

				contact.node_node_mat[contId].row[sn] = new SpVect(nnSize);

				Vect u = node.u;
				Vect v = node.getCoord().add(u);
	

				if (contact.master_edge_size[contId] == 0) {
					double length = 0;
					for (int k = 0; k < contact.master_entities[contId].length; k++) {

						int[] nids=contact.master_entities[contId][k].nodeIds;
						
						Node node1 = model.node[nids[0]];
						Node node2 = model.node[nids[1]];

						Vect v12 = node1.getCoord().add(node2.getCoord());

						length += v12.norm();
					}
					contact.master_edge_size[contId] = length;
				}

				for (int k = 0; k < contact.master_entities[contId].length; k++) {

					int[] nids=contact.master_entities[contId][k].nodeIds;
					
					Node node1 = model.node[nids[0]];
					Node node2 = model.node[nids[1]];
					int mn1 = node1.id;
					int mn2 = node2.id;
					Vect u1 = node1.u;
					Vect u2 = node2.u;
					Vect v1 = node1.getCoord().add(u1);
					Vect v2 = node2.getCoord().add(u2);
					Vect v12 = v2.sub(v1);

					v1 = v1.add(v12.times(-extention_fact));
					v2 = v2.add(v12.times(extention_fact));
					
					v12 = v2.sub(v1);
					double edgeLength = v12.norm();

					Vect edgeDir = v12.times(1. / edgeLength);

					// v12.hshow();
					Vect v1v = v.sub(v1);
					Vect v2v = v.sub(v2);

					double dot1 = v1v.dot(v12);
					double dot2 = v2v.dot(v12);
					if (dot1 != 0)
						dot1 /= v1v.norm() * edgeLength;
					if (dot2 != 0)
						dot2 /= v2v.norm() * edgeLength;
					

					
					if (dot1 * dot2 > .0) {

						contact.contacting[contId][sn] = false;

						continue;
					}

					Element elem = contact.masterElems[contId][k];
					int[] vn = elem.getVertNumb();
					for (int j = 0; j < vn.length; j++) {
						if (vn[j] != mn1 && vn[j] != mn2) {
							Node node3 = model.node[vn[j]];

							Vect v3 = node3.getCoord().add(node3.u);
							Vect v13 = v3.sub(v1);

							Vect cross1 = edgeDir.v3().cross(v13.v3());
							Vect cross2 = edgeDir.v3().cross(cross1.v3());
							contact.normals[contId][k] = cross2.normalized().v2();
							break;
						}
					}

					Vect normal = 	contact.normals[contId][k];

					double pen = v1v.dot(normal);

					if (pen < -100 * edgeLength) {
						// weakenning.el[p]=0;

						contact.contacting[contId][sn] = false;
						continue;
					}


					// if(!gradualSeperation){
					if (pen > clrFact * contact.master_edge_size[contId]) {



						contact.contacting[contId][sn] = false;
						continue;
					
					}
	

					contact.normalIndex[contId][i] = k;

					double beta = v1v.dot(edgeDir) / edgeLength;
					double alpha = 1 - beta;

					contact.node_node_mat[contId].row[sn] = new SpVect(nnSize, 2);
					contact.node_node_mat[contId].row[sn].index[0] = node1.id;
					contact.node_node_mat[contId].row[sn].index[1] = node2.id;
					contact.node_node_mat[contId].row[sn].el[0] = alpha;
					contact.node_node_mat[contId].row[sn].el[1] = beta;

	
					Vect deltaDisp = u.sub(u1.times(alpha).add(u2.times(beta)));


					//Vect deltaDisp = v.sub(v1.times(alpha).add(v2.times(beta)));
					
					Vect ur=new Vect(u.length);
					Vect u1r=new Vect(u.length);
					Vect u2r=new Vect(u.length);
						
					if(contact.frictional){
					for(int p=0;p<model.dim;p++){
					if(u_index[sn][p]>=0)
						ur.el[p]=ref_stick.el[u_index[sn][p]];
		
						int loc=u_index[mn1][p];
						if(loc>=0)
							u1r.el[p]=ref_stick.el[loc];
				
						loc=u_index[mn2][p];
						if(loc>=0)
							u2r.el[p]=ref_stick.el[loc];

					}
				}
					
					Vect vr=node.getCoord().add(ur);
					Vect v1r=node1.getCoord().add(u1r);
					Vect v2r=node2.getCoord().add(u2r);
					
				
					v1r = v1r.add(v12.times(-extention_fact));
					v2r = v2r.add(v12.times(extention_fact));
					
					Vect deltaDispRef = vr.sub(v1r.times(alpha).add(v2r.times(beta)));
					
					
				//	deltaDisp=deltaDisp.sub(deltaDispRef);
					double proj = deltaDisp.dot(normal);

					Vect disp_tang = deltaDisp.sub(normal.times(proj));
					
				//	disp_tang.show("%10.5e");

					// util.pr("-----< "+disp_tang.norm());

					Vect tang = disp_tang.normalized();
					// if(tang.norm()==0)
					tang = edgeDir.deepCopy();

					contact.tangentials[contId][k] = tang.deepCopy();
	
					contact.gap[contId].el[sn] = pen;
					contact.slide[contId].el[sn] = deltaDisp.dot(tang);

					contact.contacting[contId][sn] = true;

					totalnumContactingNodes++;

					contact.numContactingNodes[contId]++;
					break;//
				}

			}
		}

	//	resetFreedNodes();

		countStickSlip();
		/// new SpVect(weights).shownzA();
		for (int contId = 0; contId < contact.numContacts; contId++)
			util.pr((contId + 1) + ", num slave nodes: " + contact.slaveNodes[contId].length + ",  numContactingNodes: "
					+ contact.numContactingNodes[contId]);

		// new SpVect(weights).shownzA();
	}

	private void obtain_node_node3D() {

	


		totalnumContactingNodes = 0;

		int nnSize = model.numberOfNodes + 1;
		for (int contId = 0; contId < contact.numContacts; contId++) {
			
			contact.gap[contId].zero();
			contact.slide[contId].zero();	
			
			contact.numContactingNodes[contId]=0;
			

			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {
				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;

				contact.node_node_mat[contId].row[sn] = new SpVect(nnSize);

				Vect u = node.u;
				Vect v = node.getCoord().add(u);

				if (contact.master_edge_size[contId] == 0) {
					double length = 0;
					for (int k = 0; k < contact.master_entities[contId].length; k++) {

						int[] nids=contact.master_entities[contId][k].nodeIds;
						
						Node node1 = model.node[nids[0]];
						Node node2 = model.node[nids[1]];
						Node node3 = model.node[nids[2]];
						Vect v12 = node1.getCoord().add(node3.getCoord());

						double len=v12.norm();
						
						contact.master_entities[contId][k].length=len;
						length += len;
						
					}
					contact.master_edge_size[contId] = length;
				}

				for (int k = 0; k <contact. master_entities[contId].length; k++) {

					int[] nids=contact.master_entities[contId][k].nodeIds;

					int nnc=nids.length;
							
					Node[]  m_nodes= new Node[nnc];
					
					for(int m=0;m<nnc;m++)
					m_nodes[m] = model.node[nids[m]];
		
	
					Vect[] um = new Vect[nnc];
					for(int m=0;m<nnc;m++)
						um[m] = m_nodes[m].u;
				

					Vect[] vm = new Vect[nnc];
					Vect center=new Vect(3);
					for(int m=0;m<nnc;m++){
						vm[m] = m_nodes[m].getCoord().add(um[m]);
						center=center.add(vm[m]);
					}
					center.timesVoid(1./nnc);

					
					Vect v1=vm[0];
					Vect v2=vm[1];
					Vect v3=vm[2];
					Vect v4=vm[3];
		
					Vect v13=v3.sub(v1);
					Vect v24=v4.sub(v2);
					v1 = v1.add(v13.times(-extention_fact));
					v3 = v3.add(v13.times(extention_fact));
					v2 = v2.add(v24.times(-extention_fact));
					v4 = v4.add(v24.times(extention_fact));
					
					Vect cross[]=new Vect[4];
					
					cross[0]=v2.sub(v1).cross(v.sub(v1));
					cross[1]=v3.sub(v2).cross(v.sub(v2));
					cross[2]=v4.sub(v3).cross(v.sub(v3));
					cross[3]=v1.sub(v4).cross(v.sub(v4));
					
				
					
					boolean inside=true;
					for(int j=0;j<4;j++){
						for(int m=j+1;m<4;m++){
							double dot=cross[j].dot(cross[m]);
							if(dot<0) {
								inside=false;
								break;
							}
							}
						if(!inside) break;
						}
					
					if(!inside) continue;
	
					
					Vect v12=v2.sub(v1);
					Vect v23=v3.sub(v2);
					
					Vect normal = v12.cross(v23).normalized();
					
					Vect v34=v4.sub(v3);
					Vect v41=v1.sub(v4);
					
					Vect normal2 = v34.cross(v41).normalized();
			
					
					Vect normal3 = v23.cross(v34).normalized();
						
					Vect normal4 = v41.cross(v12).normalized();


					normal=normal.add(normal2).add(normal3).add(normal4);
							
					normal = normal.normalized();
				
			//	normal=new Vect(0,-1,0);
					contact.normals[contId][k] = normal.deepCopy();

				//	 normal.hshow();

					Vect v1v=v.sub(v1);
					Vect cv=v.sub(center);
					double pen = cv.dot(normal);


					if (pen > clrFact * contact.master_edge_size[contId]) {

						contact.contacting[contId][sn]= false;
						continue;
								
					}

					contact.normalIndex[contId][i] = k;

					
			
					double proj = v1v.dot(normal);

					
					double [] ww=obtainWeights( v1v, v1,  v2,  v3,  v4,  normal);
		

					 contact.node_node_mat[contId].row[sn] = new SpVect(nnSize, 4);
					for (int m = 0; m < 4; m++) {
						contact.node_node_mat[contId].row[sn].index[m] = nids[m];
						///if (m == nnn)
						contact.node_node_mat[contId].row[sn].el[m] = ww[m];

					}
					//normals[contId][k].zero();

						int p = u_index[sn][0];
					if (p < 0)
						p = u_index[sn][1];
					if (p < 0)
						p = u_index[sn][2];
					if (p < 0)
						continue;

				//	pen=u.el[2];

					contact.gap[contId].el[sn] = pen;
					
					if(contact.frictional){	

/*					Vect deltaDisp = v.deepCopy();
					for(int m=0;m<nnc;m++){
						deltaDisp=deltaDisp.sub(vm[m].times(ww[m]));
						
					}*/
					
					Vect deltaDisp = u.deepCopy();
					for(int m=0;m<nnc;m++){
						deltaDisp=deltaDisp.sub(um[m].times(ww[m]));
						
					}
					
					Vect[] umr = new Vect[nnc];	
					Vect ur = new Vect(model.dim);

					for(int q=0;q<nids.length;q++){
						umr[q]= new Vect(model.dim);
					}
					
					Vect vr;
					Vect[] vmr = new Vect[nnc];
	
		
					for(int j=0;j<model.dim;j++){
						if(u_index[sn][j]>=0)
							ur.el[j]=ref_stick.el[u_index[sn][j]];
			
						for(int q=0;q<nids.length;q++){
							
							int mn=nids[q];
							int loc=u_index[mn][j];;//util.pr(loc);util.pr(disp.length);
							if(loc>=0)
								umr[q].el[j]=ref_stick.el[loc];
					
						}

					}
					
					vr=node.getCoord().add(ur);
					
					for(int q=0;q<nids.length;q++){
						vmr[q]=m_nodes[q].getCoord().add(umr[q]);
					}
				
					
					Vect deltaDispRef = vr.deepCopy();
					for(int m=0;m<nnc;m++){
						deltaDispRef=deltaDispRef.sub(vmr[m].times(ww[m]));
						
					}
					
				//	deltaDisp=deltaDisp.sub(deltaDispRef);

					
					 proj = deltaDisp.dot(normal);
					 

					Vect disp_tang = deltaDisp.sub(normal.times(proj));
					
					double norm=disp_tang.norm();
					
					Vect tang=disp_tang.deepCopy();
					if(norm==0){
					//tang=new Vect(1,0,0).normalized();
					}
					else
						tang=disp_tang.times(1./norm);
				//	util.pr(sn);
					//tang.hshow();
				//	tang=new Vect(1,0,0).normalized();

					contact.tangentials[contId][k] = tang.deepCopy();

					//double sld=(disp_tang.dot(tang));
					double sld=Math.abs(disp_tang.dot(tang));

					contact.slide[contId].el[sn]=sld;
					}
					
					contact.contacting[contId][sn] = true;

					totalnumContactingNodes++;

					contact.numContactingNodes[contId]++;
					break;//
				}

			}
		}

	 //node_node_mat[contId].shownzA();

		//resetFreedNodes();

		//countStickSlip();
		/// new SpVect(weights).shownzA();
		for (int contId = 0; contId < contact.numContacts; contId++)
			util.pr((contId + 1) + ", num slave nodes: " + contact.slaveNodes[contId].length + ",  numContactingNodes: "
					+ contact.numContactingNodes[contId]);

		//
	///new SpVect(gap).shownzA();
	}

	private void updateGap(boolean allowSep) {


		for (int contId = 0; contId < contact.numContacts; contId++) {

			contact.gap[contId].zero();

			contact.slide[contId].zero();

			
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {
				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;

				if (!contact.contacting[contId][sn])
					continue;

				Vect u = node.u;
				Vect v = node.getCoord().add(u);

				if (contact.node_node_mat[contId].row[sn].index == null)
					continue;

				int mn1 = contact.node_node_mat[contId].row[sn].index[0];
				int mn2 = contact.node_node_mat[contId].row[sn].index[1];

				Node node1 = model.node[mn1];
				Node node2 = model.node[mn2];
				Vect u1 = node1.u;
				Vect u2 = node2.u;
				Vect v1 = node1.getCoord().add(u1);
				Vect v2 = node2.getCoord().add(u2);
				Vect v12 = v2.sub(v1);

				v1 = v1.add(v12.times(-extention_fact));
				v2 = v2.add(v12.times(extention_fact));

				v12 = v2.sub(v1);

				Vect v1v = v.sub(v1);

				Vect edgeDir = v2.sub(v1).normalized();
				double edgeLength = v12.norm();

				int nrmIndex = contact.normalIndex[contId][i];

				Vect normal = null;
				Element elem = contact.masterElems[contId][nrmIndex];
				int[] vn = elem.getVertNumb();
				for (int j = 0; j < vn.length; j++) {
					if (vn[j] != mn1 && vn[j] != mn2) {
						Node node3 = model.node[vn[j]];

						Vect v3 = node3.getCoord().add(node3.u);
						Vect v13 = v3.sub(v1);

						Vect cross1 = edgeDir.v3().cross(v13.v3());
						Vect cross2 = edgeDir.v3().cross(cross1.v3());
						normal = cross2.normalized().v2();
						break;
					}
				}

				double pen = v1v.dot(normal);

				if (allowSep && pen > clrFact * contact.master_edge_size[contId]) {

					continue;
				}

				if (pen < -100 * edgeLength) {

					continue;
				}

				double alpha = contact.node_node_mat[contId].row[sn].el[0];
				double beta = contact.node_node_mat[contId].row[sn].el[1];

				Vect deltaDisp = u.sub(u1.times(alpha).add(u2.times(beta)));
				/// deltaDisp.times(1e9).hshow();
				double proj = deltaDisp.dot(normal);

				Vect disp_tang = deltaDisp.sub(normal.times(proj));

				Vect tang = disp_tang.normalized();

				tang = edgeDir.deepCopy();

				int p = u_index[sn][0];
				if (p < 0)
					p = u_index[sn][1];
				if (p < 0)
					continue;

				contact.gap[contId].el[sn] = pen;

				contact.slide[contId].el[sn] = deltaDisp.dot(tang);

			}
		}


	}

	private void assembleConstraintMats() {

		int dof = model.Ks.nRow;
		int nRows=model.numberOfNodes + 1;
		
		for (int contId = 0; contId < contact.numContacts; contId++){
		contact.constraint_matrix_N[contId] = new SpMatAsym(nRows, dof);
		contact.constraint_matrix_T[contId] = new SpMatAsym(nRows, dof);

		contact.constraint_matrix_STK[contId] = new SpMatAsym(nRows, dof);
		}


		for (int contId = 0; contId < contact.numContacts; contId++){
			
		
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {

				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;
	

				if (contact.node_node_mat[contId].row[sn].nzLength > 0) {

					int index = model.U_unknownIndex[sn] - 1;
					if (index < 0)
						continue;

					int px = u_index[sn][0];
					int py = u_index[sn][1];
					int pz = -1;
					if (model.dim == 3)
						pz = u_index[sn][2];

	

					contact.constraint_matrix_N[contId].row[sn] = new SpVect(dof);

					Vect normal = contact.normals[contId][contact.normalIndex[contId][i]];

				
					int[] nids=contact.master_entities[contId][contact.normalIndex[contId][i]].nodeIds;

						
					int mn1 = nids[0];
					int mn2 = nids[1];

					int p1x = u_index[mn1][0];
					int p1y = u_index[mn1][1];

					int p2x = u_index[mn2][0];
					int p2y = u_index[mn2][1];

					double alpha = contact.node_node_mat[contId].row[sn].el[0];
					double beta = contact.node_node_mat[contId].row[sn].el[1];

					contact.constraint_matrix_N[contId].row[sn] = new SpVect(dof, 6);

					int kx = 0;
					if (px != -1) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = px;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = normal.el[0];
					}
					if (py != -1) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = py;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = normal.el[1];
					}

					// util.pr(weights.el[com_index]);

					if (p1x >= 0) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = p1x;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = -alpha * normal.el[0];
					}
					if (p1y >= 0) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = p1y;
						;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = -alpha * normal.el[1];
					}

					if (p2x >= 0) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = p2x;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = -beta * normal.el[0];
					}
					if (p2y >= 0) {
						contact.constraint_matrix_N[contId].row[sn].index[kx] = p2y;
						;
						contact.constraint_matrix_N[contId].row[sn].el[kx++] = -beta * normal.el[1];
					}
					contact.constraint_matrix_N[contId].row[sn].sortAndTrim(kx);
					

					// Vect tang=new Vect(-normal.el[1],normal.el[0]);

					Vect tang = contact.tangentials[contId][contact.normalIndex[contId][i]].deepCopy();

				////	tang.hshow();
				
					// ===

					if (contact.fric_coef[contId] == 0) {
						tang.zero();
					}
					contact.constraint_matrix_T[contId].row[sn] = new SpVect(dof, 6);
					kx = 0;
					if (px != -1) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = px;
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = tang.el[0];
					}
					if (p1y != -1) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = py;
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = tang.el[1];
					}

					if (p1x >= 0) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = p1x;
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = -alpha * tang.el[0];
					}
					if (p1y >= 0) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = p1y;
						;
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = -alpha * tang.el[1];
					}

					if (p2x >= 0) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = p2x;
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = -beta * tang.el[0];
					}
					if (p2y >= 0) {
						contact.constraint_matrix_T[contId].row[sn].index[kx] = p2y;
						
						contact.constraint_matrix_T[contId].row[sn].el[kx++] = -beta * tang.el[1];
					}

					if (contact.stick[contId][sn]) {
						contact.constraint_matrix_STK[contId].row[sn] = contact.constraint_matrix_T[contId].row[sn].deepCopy();

					} else {
						contact.constraint_matrix_STK[contId].row[sn] = contact.constraint_matrix_T[contId].row[sn].times(reduct);
					}

				}
			}
				
			//	constraint_matrix_TX.addGeneral(contact.constraint_matrix_T[contId]);
			//	contact.constraint_matrix_T[contId].shownzA();
			}
		
	//	util.pr("GcfTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "+constraint_matrix_TX.norm2());


		//constraint_matrix_TX.shownzA();
		
		//util.pr("constraint_matrix_T_ norm2 "+constraint_matrix_TX.norm2());


		// G_stk.shownzA();

	}

	private void assembleConstraintMats3D() {

		int dof = model.Ks.nRow;
		int nRows=model.numberOfNodes + 1;
		
		for (int contId = 0; contId < contact.numContacts; contId++){
		contact.constraint_matrix_N[contId] = new SpMatAsym(nRows, dof);
		contact.constraint_matrix_T[contId] = new SpMatAsym(nRows, dof);

		contact.constraint_matrix_STK[contId] = new SpMatAsym(dof, dof);
		}

		for (int contId = 0; contId <  contact.numContacts; contId++)
			for (int i = 0; i <  contact.slaveNodes[contId].length; i++) {

				Node node =  contact.slaveNodes[contId][i];
				int sn = node.id;
				

				if ( contact.node_node_mat[contId].row[sn].nzLength > 0) {

					int index = model.U_unknownIndex[sn] - 1;
					if (index < 0)
						continue;
					

						int px = u_index[sn][0];
						int py = u_index[sn][1];
						int pz = u_index[sn][2];

						int p = px;
						if (p == -1)
							p = py;
						if (p == -1)
							p = pz;


						Vect normal = contact.normals[contId][ contact.normalIndex[contId][i]];

					
						int[] nids =  contact.master_entities[contId][ contact.normalIndex[contId][i]].nodeIds;
								
					
						int ps[]=new int[3];

						for(int k=0;k<3;k++)
						ps[k] = u_index[sn][k];
						
						
						int pm[][]=new int[nids.length][3];
						for(int j=0;j<nids.length;j++)
						for(int k=0;k<3;k++)
							pm[j][k] = u_index[nids[j]][k];


						contact.constraint_matrix_N[contId].row[sn] = new SpVect(dof, 15);
						
						double wgt[] =  contact.node_node_mat[contId].row[sn].el;


						int kx = 0;
						
						for(int k=0;k<3;k++)
							if(ps[k]!=-1){
								contact.constraint_matrix_N[contId].row[sn].index[kx] = ps[k];
								contact.constraint_matrix_N[contId].row[sn].el[kx++] = normal.el[k];
							}
		
						// util.pr(weights.el[com_index]);
						for(int j=0;j<nids.length;j++){
							for(int k=0;k<3;k++)
								if(pm[j][k]!=-1){
								contact.constraint_matrix_N[contId].row[sn].index[kx] = pm[j][k];
							contact.constraint_matrix_N[contId].row[sn].el[kx++] = -wgt[j] * normal.el[k];
								}
								
						}

						contact.constraint_matrix_N[contId].row[sn].sortAndTrim(kx);

						// Vect tang=new Vect(-normal.el[1],normal.el[0]);

						if(contact.frictional){
						Vect tang = contact.tangentials[contId][ contact.normalIndex[contId][i]].deepCopy();

					

						if ( contact.fric_coef[contId] == 0) {
							tang.zero();
						}
						contact.constraint_matrix_T[contId].row[sn] = new SpVect(dof, 15);
						
						//==================================
						kx=0;
						for(int k=0;k<3;k++){
							if(ps[k]!=-1){
								contact.constraint_matrix_T[contId].row[sn].index[kx] = ps[k];
								contact.constraint_matrix_T[contId].row[sn].el[kx++] = tang.el[k];
							
		
						// util.pr(weights.el[com_index]);
						for(int j=0;j<nids.length;j++){
								if(pm[j][k]!=-1){
								contact.constraint_matrix_T[contId].row[sn].index[kx] = pm[j][k];
								contact.constraint_matrix_T[contId].row[sn].el[kx++] = -wgt[j] * tang.el[k];
								}
							}
						
						}
							
						}
						
						contact.constraint_matrix_T[contId].row[sn].sortAndTrim(kx);
		
						contact.constraint_matrix_STK[contId].row[sn] = new SpVect(dof, 15);

						
						int[] kk=new int[3];
						
						for(int k=0;k<3;k++){
							if(ps[k]!=-1 && contact.stick[contId][sn]){
								contact.constraint_matrix_STK[contId].row[ps[k]] = new SpVect(dof, 10);

								contact.constraint_matrix_STK[contId].row[ps[k]].index[kk[k]] = ps[k];
								contact.constraint_matrix_STK[contId].row[ps[k]].el[kk[k]++] = 1-normal.el[k];
							
	
						for(int j=0;j<nids.length;j++){
								if(pm[j][k]!=-1){
								contact.constraint_matrix_STK[contId].row[ps[k]].index[kk[k]] = pm[j][k];
							contact.constraint_matrix_STK[contId].row[ps[k]].el[kk[k]++] = -wgt[j] * (1-normal.el[k]);
								}
								
						}
						
						contact.constraint_matrix_STK[contId].row[ps[k]].sortAndTrim(kk[k]);
					}
						}
						}
					
					
				
			}
		}
		//util.pr(" here 555");
		// G_stk.shownzA();
	}

	private void checkStickSlip() {

		for (int contId = 0; contId < contact.numContacts; contId++) {
			double mu = contact.fric_coef[contId];
			if (mu == 0)
				continue;
			
			
			
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {

				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;

				if (!contact.contacting[contId][sn]) {
					contact.stick[contId][sn] = false;
					contact.landed_stick[contId][sn] = false;
					continue;
				}

				int index = model.U_unknownIndex[sn] - 1;
				if (index < 0)
					continue;
				
				if (contact.lamN[contId].el[sn] >=0)  {
					contact.stick[contId][sn] = false;
					contact.landed_stick[contId][sn] = false;
					contact.lamT[contId].el[sn] = 0;
					contact.slide[contId].el[sn] = 0;
					contact.contacting[contId][sn] = false;
					continue;
				}

				if (contact.node_node_mat[contId].row[sn].nzLength > 0) {
					if (contact.lamN[contId].el[sn] < 0) {
						double abs_lamT = Math.abs(contact.lamT[contId].el[sn]);
						double muFn = mu * Math.abs(contact.lamN[contId].el[sn]);
			///	util.pr(contact.lamT[contId].el[sn] +"   "+contact.lamN[contId].el[sn]);
						if (abs_lamT > muFn * (1 + margin)) {
							if (contact.lamT[contId].el[sn] > 0)
								contact.lamT[contId].el[sn] = muFn;
							else
								contact.lamT[contId].el[sn] = -muFn;
							contact.stick[contId][sn] = false;
							contact.landed_stick[contId][sn] = false;
							contact.slide[contId].el[sn] = 0;
						} else{
							if (!contact.stick[contId][sn]) {
								contact.stick[contId][sn] = true;

								contact.landed_stick[contId][sn] = true;
								contact.slide[contId].el[sn] = 0;
							}
						}
					} 	if (contact.lamN[contId].el[sn] >0)  {
						contact.stick[contId][sn] = false;
						contact.landed_stick[contId][sn] = false;
						contact.lamT[contId].el[sn] = 0;
						contact.slide[contId].el[sn] = 0;
						contact.contacting[contId][sn] = false;
					}
				} else {
					contact.stick[contId][sn] = false;
					contact.landed_stick[contId][sn] = false;
					contact.lamT[contId].el[sn] = 0;
					contact.slide[contId].el[sn] = 0;
					contact.contacting[contId][sn] = false;
				}
			}

		}
		
		for (int contId = 0; contId < contact.numContacts; contId++) {
			double mu =  contact.fric_coef[contId];
			if (mu == 0)
				continue;
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {

				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;


			if (contact.landed_stick[contId][sn]) {
				
				int[] nids=contact.master_entities[contId][contact.normalIndex[contId][i]].nodeIds;


				for(int p=0;p<model.dim;p++){
				if(u_index[sn][p]>=0)
					ref_stick.el[u_index[sn][p]]=disp.el[u_index[sn][p]];
				
				for(int q=0;q<nids.length;q++){
					int mn=nids[q];
					int loc=u_index[mn][p];
					if(loc>=0)
						ref_stick.el[loc]=disp.el[loc];
				}
				}
		
				contact.landed_stick[contId][sn] = false;

			}
		}
		
			}
		
		countStickSlip();

		//new SpVect(ref_stick).shownzA();
		
		 
		/// G_stkt.shownzA();

		// G_stk.shownzA();

	}

	private void printoutErrors(){

		util.pr("NR error");

		int nn = 0;
		for (int k = 0; k < nr_err.length; k++)
			if (nr_err.el[k] > 0)
				nn++;

		Vect nr_distilled = new Vect(nn);

		int nr_it_distilled[] = new int[nn];
		int kx = 0;
		for (int k = 0; k < nr_err.length; k++) {
			if (nr_err.el[k] > 0) {
				nr_it_distilled[kx] = nr_it[k];

				nr_distilled.el[kx++] = (nr_err.el[k]);
			}
		}

		// nr_distilled.show("%5.4e");
		for (int k = 0; k < nr_distilled.length; k++) {
			System.out.format(" %d\t%6.4e\n", nr_it_distilled[k], nr_distilled.el[k]);
			nr_distilled.el[k] = Math.log10(nr_distilled.el[k]);
		}
	
		util.plot(nr_distilled);
	

		util.pr("aug_disp changes vs aug_iter");
		aug_disp_err.show("%5.4e");

		util.pr("Gap[m] vs aug_iter");
		gap_err.show("%5.4e");
		// util.plot(err);

		if(contact.frictional){
		util.pr("slide[m] vs aug_iter");
		slide_err.show("%5.4e");
		}
	
	}
	
	private Vect solveLinear(SpMatSolver solver, SpMat Ks, Vect dF) {

		Vect du = new Vect(dF.length);

		model.Ci = Ks.scale(dF);

		// if(model.Ls==null)
		model.Ls = Ks.ichol();

		if (dF.abs().max() != 0) {
			 if(model.xp==null){
			du = solver.ICCG(Ks, model.Ls, dF, model.errCGmax, model.iterMax);
			 }
			 else{
			// du=solver.ICCG(Ks,model.Ls,
			// dF,model.errCGmax,model.iterMax,model.xp);
			 du=model.solver.err0ICCG(Ks,model.Ls, dF,model.errCGmax*1e-3,model.iterMax,model.xp);

			 }
		} else {
			util.pr("Solution is zero!");
			du = new Vect(Ks.nRow);
		}

		// model.xp=du.deepCopy();

		du.timesVoid(model.Ci);

		return du;
	}

	private Vect calcResidual(SpMat Ks, Vect u, Vect b, int iter) {

		Vect Fint = K_hat.smul(u);

		Vect dF = b.sub(Fint);

		
		Fc.zero();
		Fcf.zero();
		

		for (int contId = 0; contId < contact.numContacts; contId++) {
		Vect pg = contact.gap[contId].deepCopy();

		for (int k = 0; k < pg.length; k++) {
			pg.el[k] *= contact.weights[contId].el[k];
					
		}
		
		

		if(contact.constraint_matrix_N_trp!=null)
		Fc = Fc.add(contact.constraint_matrix_N_trp[contId].mul(pg));



	    if(contact.constraint_matrix_STK_trp!=null){
	    	
	    }

	   		Vect ps=contact.constraint_matrix_STK[contId].mul(disp);

	   		ps.timesVoid(fr);
	

		for (int k = 0; k < ps.length; k++) {
			ps.el[k] *=contact.weights[contId].el[u_index_inv[k]];

		}
		
	//	new SpVect(ps).shownzA();

	
				Fcf = Fcf.add(contact.constraint_matrix_STK_trp[contId].mul(ps));
		
			

		}

		
		dF = dF.sub(Fc).sub(Fcf); //

		dF = dF.sub(aug_N).sub(aug_T);
		
		///util.pr(" aug_TTTTTTTTTTTTTTTTTTTTTT   "+aug_T.norm());


		return dF;

	}

	public void readContact(Loader loader, BufferedReader br, Model model) throws IOException {

		model.setEdge();
		

		double minEdgeLenghth = model.minEdgeLength;
		util.pr(" minEdgeLength ------------------------- " + minEdgeLenghth);
		double clearFact = 1e-4;
		
		String line = loader.getNextDataLine(br," / * N_CONTACT * NR_ITR. * AUG_ITER  * NLOADS  * MNR /");

		int[] integs=loader.getCSInt(line);
	
		int numCont = integs[0];

		nr_itmax = nr_itmax0;
		aug_itmax =aug_itmax0;
		nLoads = nLoads0;
		n_modifNR = n_modifNR0;
		if(integs.length>1){
			nr_itmax = integs[1];
		}
		
		if(integs.length>2){
			aug_itmax = integs[2];
		}
		
		if(integs.length>3){
			nLoads = integs[3];
		}
		
		if(integs.length>4){
			n_modifNR = integs[4];
		}

		
		contact=new Contact(numCont);

	
		
		nr_tol=nr_tol0;
		aug_tol=aug_tol0;
		
		mnr_tol=0;
		if(integs.length>1){
			
			 line = loader.getNextDataLine(br," / * NR_TOL  *  AUG_TOL * MNR_TOL/");
			 double[] tols=loader.getTabedData(line);
			 nr_tol=tols[0];
			 if(integs.length>2)
			 aug_tol=tols[1];
			 if(tols.length>2)
				 mnr_tol=tols[2];
		}
		
		gap_tol=aug_tol;
		
		if(mnr_tol==0) mnr_tol=nr_tol;
		
		contact.readContacts(loader, br, model);

		

	}

	


	private void calcPenaltyFactor() {

		penalMax = 0;

		for (int contId = 0; contId < contact.numContacts; contId++)
			for (int i = 0; i < contact.slaveNodes[contId].length; i++) {
				int sn = contact.slaveNodes[contId][i].id;

				int index = model.U_unknownIndex[sn] - 1;

				if (index < 0)
					continue;
				int xind = u_index[sn][0];
				int yind = u_index[sn][1];
				int zind = -1;
				if (model.dim == 3)
					zind = u_index[sn][2];
				double val1 = 0;
				double val2 = 0;
				double val3 = 0;

				if (xind != -1)
					val1 =  model.Ks.row[xind].el[model.Ks.row[xind].nzLength - 1];

				if (yind != -1)
					val2 =  model.Ks.row[yind].el[model.Ks.row[yind].nzLength - 1];
				if (zind != -1)
					val3 =  model.Ks.row[zind].el[model.Ks.row[zind].nzLength - 1];


				// double max=(val1+val2)/2;
				double max = (val1+val2+val3);
				
	
			/////	if(model.dim==3) max/=10;

				if (max > penalMax)
					penalMax = max;

					contact.weights[contId].el[sn] = 	contact.penFactor[contId] * max;


				// weightsf.el[p]=weights.el[p]*fn_ratio[contId];

				////// weights.el[p]*=1./slaveNodes[contId].length;

			}

		// weights0=weights.deepCopy();



		util.pr("pf max:" + penalMax);
		util.pr("pft max:" + fr*penalMax);


		
	}

	private void countStickSlip() {

		for (int contId = 0; contId < 	contact.slaveNodes.length; contId++) {
			int nstk = 0;
			int nslip = 0;
			for (int i = 0; i < 	contact.slaveNodes[contId].length; i++) {
				Node node = 	contact.slaveNodes[contId][i];
				int sn = node.id;

				if (contact.contacting[contId][sn]) {


					if (contact.stick[contId][sn])
						nstk++;
					else
						nslip++;

				}
				// else
				// util.pr("node "+sn+" free ");
			}
			util.pr((contId + 1) + ", stick: " + nstk + ",  slip: " + nslip);
		}
	}

	private boolean handleModifiedMNR(int nr_iter, SpMatSolver solver, boolean direct)
	{
		
		double er=1;
		boolean modif_done = false;

		if ( /*disp_err<nr_tol &&  */totalNRIter > 0 && n_modifNR > 0) {
			if (direct && n_modifNR > 0) {
				KK = Ks.matForm();
				KK.lu();
			}

			er=modifiedNR(solver, rhs, direct, nr_tol);

			modif_done = true;
		}

		if (modif_done) {
			
			if (mnr_tol<nr_tol && er < nr_tol) {
				nr_err.el[totalNRIter] = er;
				nr_it[totalNRIter] = totalNRIter;
				util.pr("nr_iter: " + nr_iter + "           nr_err: " + er + "       disp_err: ");

				return true;
			}
			assembleContactMatrices();
			addMatrices();

		}
		
		return false;
	}
	
	private void updateLambdas(){
		for(int contId=0;contId<contact.numContacts;contId++){
			
			Vect gap=contact.gap[contId];
			Vect slide=contact.slide[contId];

			Vect weights=contact.weights[contId];


		double relax=0.5; // it should be 1.0 ideally

		Vect pgap = gap.deepCopy();
		Vect pslide = slide.deepCopy();

		for (int k = 0; k < gap.length; k++) {
			pgap.el[k] *= weights.el[k];
		
			pslide.el[k] *= weights.el[k] * fr;
		}


		for (int count = 0; count < contact.slaveNodes[contId].length; count++) {
			
			int snid=contact.slaveNodes[contId][count].id;


			if(!contact.contacting[contId][snid]) continue;
			
			
			if(gap.el[snid]>=0){
				gap.el[snid]=0;
				contact.lamN[contId].el[snid] = 0;
				contact.lamT[contId].el[snid] = 0;
				contact.contacting[contId][snid]=false;
				contact.stick[contId][snid]=false;
				continue;
			}
			

			
				double aug =contact.lamN[contId].el[snid] +pgap.el[snid]*relax;

				contact.lamN[contId].el[snid] = aug ;
			
				double augT =contact.lamT[contId].el[snid] +pslide.el[snid]*relax;

				contact.lamT[contId].el[snid] = augT;


		}

	
		}
	}
	
	private void updateAugLoads(){

		aug_N.zero();
		aug_T.zero();
		


		for(int contId=0;contId<contact.numContacts;contId++){
			
			
			Vect lamN=contact.lamN[contId];
			Vect lamT=contact.lamT[contId];
			
		if (!aug_normal)
			lamN.zero();

		if (!aug_tang)
			lamT.zero();

		if(contact.constraint_matrix_N_trp[contId]!=null)
		aug_N = aug_N.add(contact.constraint_matrix_N_trp[contId].mul(lamN));
		
		if(contact.constraint_matrix_T_trp[contId]!=null){
		aug_T =  aug_T.add(contact.constraint_matrix_T_trp[contId].mul(lamT));

		}
		}
	}
	
	private void writeContactForce(){
		


		for (int contId = 0; contId < contact.numContacts; contId++)
			for (int k = 0; k < contact.master_entities[contId].length; k++) {

			
				
				int[] nids = contact.master_entities[contId][k].nodeIds;
	

				Vect normal = contact.normals[contId][k];
				
				if (normal == null)
					continue;
				
			//	Vect normal1=normal.deepCopy();

				//normal1.el[1]=0;
				
				for(int m=0;m<nids.length;m++){
				model.node[nids[m]].Fms = normal.deepCopy();
	
				for(int p=0;p<model.dim;p++){
					model.node[nids[m]].Fms.el[p]+=1e-6*Math.random();
				
				}
				}
				// node2.u=normal.deepCopy();
			}

	model.writeNodalField(model.resultFolder + "\\results\\master_normal.txt", 2);

	for (int i = 1; i <= model.numberOfNodes; i++) {
		Vect F = model.node[i].Fms;
		if (F != null)
			model.node[i].Fms = null;
	}

	for (int contId = 0; contId < contact.numContacts; contId++)
		for (int k = 0; k < contact.slaveNodes[contId].length; k++) {

			Node node = contact.slaveNodes[contId][k];

			int ind = contact.normalIndex[contId][k];

			if (ind < 0)
				continue;

			Vect normal = contact.normals[contId][ind];// new
												// Vect(-v12.el[1],v12.el[0]).normalized();
			if (normal == null)
				continue;

			model.node[node.id].Fms = normal.times(-1);
			for(int p=0;p<model.dim;p++){
				model.node[node.id].Fms.el[p]+=1e-6*Math.random();
			}
			// node2.u=normal.deepCopy();
		}
	model.writeNodalField(model.resultFolder + "\\results\\slave_normal.txt", 2);

	Fcf=Kcf.smul(disp);
///	new SpVect(Fc).shownzA();
	if(aug_N.norm()==0)
		model.setU(Fc.times(1).add(Fcf.times(1)).times(-1));
	else
	 model.setU(aug_N.times(1).add(aug_T.times(1)).times(-1));
	for (int i = 1; i <= model.numberOfNodes; i++) {
		// model.node[i].u.el[1]=0;
	}
	
	model.writeNodalField(model.resultFolder + "\\results\\contact_force.txt", -1);
	
	}
	
	
	private void plotRadialDisp(Vect xr, Vect[] urs, int[] xr_nids, int aug_iter){
		for (int k = 0; k < xr.length; k++) {
			// int n=2551+k;
			int n = xr_nids[k];

			int index = model.U_unknownIndex[n] - 1;
			if (index < 0)
				continue;

			int com_index = u_index[n][0];
			// urs[aug_iter].el[k]=-Math.abs(u.el[com_index])*1e6;
			if (com_index == -1)
				continue;
			urs[aug_iter].el[k] = disp.el[com_index] * 1e3;
			if (xr.el[k] < 0)
				urs[aug_iter].el[k] *= -1;

			//util.pr(urs[aug_iter].el[k]);
		}
	}
	
	private void initialize() {
		
		
		contact.initialize(model);
		

		int[][] ne = new int[model.numberOfNodes + 1][20];
		int[] nz = new int[model.numberOfNodes + 1];

		for (int i = 1; i <= model.numberOfElements; i++) {
			int[] vn = model.element[i].getVertNumb();
			for (int j = 0; j < vn.length; j++) {
				ne[vn[j]][nz[vn[j]]++] = i;

			}
		}


		u_index = new int[model.numberOfNodes + 1][model.dim];
		u_index_inv=new int[model.numberOfNodes*model.dim];
		for (int i = 1; i <= model.numberOfNodes; i++)
			for (int k = 0; k < model.dim; k++)
				u_index[i][k] = -1;

		int ix = 0;
		for (int i = 1; i <= model.numberOfUnknownU; i++) {
			int nodeNumb = model.unknownUnumber[i];

			if (model.node[nodeNumb].isDeformable()) {
				for (int k = 0; k < model.dim; k++) {
					if (!model.node[nodeNumb].is_U_known(k)) {
						u_index[nodeNumb][k] = ix;
						u_index_inv[ix]=i;
						ix++;
					}
				}

			}
		}


		ref_stick = new Vect(model.Ks.nRow);

		aug_N = new Vect(model.Ks.nRow);
		aug_T = new Vect(model.Ks.nRow);



		Fc = new Vect(model.Ks.nRow);
		Fcf = new Vect(model.Ks.nRow);

		top=new Vect(model.nTsteps);
		
		 gap_err = new Vect(aug_itmax*model.nTsteps);
		 aug_disp_err = new Vect(aug_itmax*model.nTsteps);
		 slide_err = new Vect(aug_itmax*model.nTsteps);

		 nr_err = new Vect(nLoads * aug_itmax * (nr_itmax + n_modifNR)*model.nTsteps);
		 nr_it = new int[nLoads * aug_itmax * (nr_itmax + n_modifNR)*model.nTsteps];
	}
	
		
	
	public Vect getDeformation(Model model, SpMatSolver solver, int mode,int step) {

		solver.terminate(false);
		System.out.println(" Static analysis....");

		double loadFactor0 = 1;

		if (model.dim == 3)
			loadFactor0 = 1;
		
		if (model.centrigForce) {
			model.setNodalMass();
			double rps = model.rpm / 30 * PI;
			double omeag2 = rps * rps;
			for (int i = 1; i <= model.numberOfNodes; i++) {
				Vect v = model.node[i].getCoord();
				double m = model.node[i].getNodalMass();
				Vect F = v.times(omeag2 * m/loadFactor0);
				model.node[i].setF(F);
			}

		}
		if(model.pressLoads!=null){
			for(int i=0;i<model.pressLoads.length;i++)
				model.pressLoads[i].setPressure(model);
		}


		Vect bU1=model.bU.add(model.getbUt(mode));
		
		bU1.timesVoid((step+1.)/model.nTsteps);
		


		boolean axi = (model.struc2D == 2);

		if (axi)
			loadFactor0 *= 2 * PI;
		
		bU1.timesVoid(loadFactor0);


////new SpVect(bU1).shownzA();
		Vect u=solve(model, solver,model.Ks,bU1,step);
		
		return u;
	}
	
	public Vect getVibration(Model model, SpMatSolver solver, int mode,int step) {
		
		
		solver.terminate(false);
		System.out.println(" Calculating vibration using the Newmark method....");

		double dt=model.dt;
		
		double loadFactor0 = 1000;
		
		if (model.dim == 3)
			loadFactor0 = 1;

		boolean axi = (model.struc2D == 2);

		if (axi)
			loadFactor0 *= 2 * PI;
		
		
		if (model.centrigForce) {
			model.setNodalMass();
			double rpm = 7000;
			double rps = rpm / 30 * PI;
			double omeag2 = rps * rps;
			for (int i = 1; i <= model.numberOfNodes; i++) {
				Vect v = model.node[i].getCoord();
				double m = model.node[i].getNodalMass();
				Vect F = v.times(omeag2 * m);
				model.node[i].setF(F);
			}

		}

		if(step==0){
			
			for (int i = 1; i <= model.numberOfNodes; i++) {
			
				Vect F = model.node[i].F;
				if(F!=null)
				model.node[i].setF(F.times(loadFactor0));
			}
			


		}

		Vect bU1=model.bU.add(model.getbUt(mode));
		
		if(step<2){
			

			
			Vect u=solve( model, solver,model.Ks ,bU1,step);
			
			if(step==1)
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
		
		///if(step==0) model.up=new Vect(bU1.length);
		
		///if(step>9) bU1.zero();

		
		SpMat Ks=model.Ks.addNew(model.Ms.timesNew(b1)).addNew(model.Cs.timesNew(b4));

			Vect bp=model.Ms.smul(model.up.times(b1).add(model.ud.times(-b2)).add(model.udd.times(-b3)))
			.add(model.Cs.smul(model.up.times(b4).add(model.ud.times(-b5)).add(model.udd.times(-b6))));

		
			bU1=bU1.add(bp);

			Vect u=solve( model, solver,Ks ,bU1,step);



		model.up=u.deepCopy();
	


		return u;
		

	}

	private double[] getGapAndSlideErrors(){
		double [] errors=new double [2];
		
		double e1max=0;
		double e2max=0;
	
		for (int contId = 0; contId < contact.numContacts; contId++) {
			
			double master_size_sq=Math.pow(contact.master_edge_size[contId],2);
			

			double sum_er1_sq=0;
			double sum_er2_sq=0;
			for (int i = 0; i <contact.slaveNodes[contId].length; i++) {
				Node node = contact.slaveNodes[contId][i];
				int sn = node.id;
				
				double pen=contact.gap[contId].el[sn];
				double sld=contact.slide[contId].el[sn];
			//	util.pr(sld);
				
				sum_er1_sq+=pen*pen;
				sum_er2_sq+=sld*sld;
			}
				double rel_er1=Math.sqrt(sum_er1_sq/master_size_sq);
				double rel_er2=Math.sqrt(sum_er2_sq/master_size_sq);
				
				if(rel_er1>e1max) e1max=rel_er1;
				if(rel_er2>e1max) e2max=rel_er2;
				
				
			}
		
		errors[0]=e1max;
		errors[1]=e2max;
		
		
	
		return errors;
	}
	
	private double [] obtainWeights(Vect v1v,Vect v1, Vect v2, Vect v3, Vect v4, Vect normal){
		

		
		double proj = v1v.dot(normal);

		Vect v1v_proj = v1v.sub(normal.times(proj));
		
		Mat R=util.rotMat(new Vect(0,0,1), normal);
		v1=R.mul(v1);
		v2=R.mul(v2);
		v3=R.mul(v3);
		v4=R.mul(v4);
		
		v1v_proj=R.mul(v1v_proj);
		
		
		Vect sn_proj = v1.add(v1v_proj);


		double a0=v2.el[0]-v1.el[0];
		double a1=v4.el[0]-v1.el[0];
		double a2=v1.el[0]-v2.el[0]+v3.el[0]-v4.el[0];

		double b0=v2.el[1]-v1.el[1];
		double b1=v4.el[1]-v1.el[1];
		double b2=v1.el[1]-v2.el[1]+v3.el[1]-v4.el[1];

		Vect P=sn_proj.sub(v1).v2();
		
		
	///	P=v1.add(v3).times(.5).v2();
		
		Vect x=new Vect(2);
		
		if(P.norm()!=0){
			
							
		Mat M1=new Mat(2,2);
		Mat M2=new Mat(2,2);


		Vect dx=new Vect(2);
		double err=1;
		for(int j=0;j<10;j++){

			if(err<1e-4) break;
		M1.el[0][0]=a0;
		M1.el[0][1]=a1+a2*x.el[0];
		M1.el[1][0]=b0;
		M1.el[1][1]=b1+b2*x.el[0];
		
	
		
		Vect b=P.sub(M1.mul(x));
		

		err=b.norm()/P.norm();

		M2.el[0][1]=a2*x.el[0];
		M2.el[1][1]=b2*x.el[0];
		
		Mat M=M1.add(M2);
		
	//	M.show();
		
		Mat invM=M.inv2();
	
	
		dx=invM.mul(b);
		x=x.add(dx);
	//	x.show();

	//	util.pr(" err  ========>   "+err);
		}
		}
		
		 double uu=x.el[0];
		 double vv=x.el[1];
		 
		 double[] ww=new double[4];
	
		 ww[0]=(1-uu)*(1-vv);
		 ww[1]=(uu)*(1-vv);
		 ww[2]=(uu*vv);
		 ww[3]=(1-uu)*(vv);
		 
		 return ww;
	}
	
	public static void main(String[] args) {

		new Main();
	}

}

