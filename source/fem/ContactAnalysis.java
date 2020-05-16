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

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan. Created Aug 20, 2012.
 */

public class ContactAnalysis {
	
	private class MasterEntity{
		
		int[] nodeIds;
		double length;
		
		private MasterEntity(int n){
			nodeIds=new int[n];
		}
	}

	private SpMat Ks = null;

	private SpMatAsym node_node = null;
	private SpMatAsym Gc = null;
	private SpMatAsym Gct = null;
	private SpMatAsym Gcf = null;
	private SpMatAsym Gcfadh = null;
	private SpMatAsym Gcft = null;
	private SpMatAsym G_stk = null;
	private SpMatAsym G_stkt = null;

	private SpMat Kc = null;
	private SpMat Kcf = null;
	private SpMat Kadh = null;
	private SpMat Kadhf = null;

	private Vect lamN, gap;
	private Vect lamT, slide;// ,slide_prev;

	Mat KK = null;
	MatSolver direct_slv = null;

	private Vect aug_N;
	private Vect aug_T;

	private Vect ref_stick;

	private Vect Fc;
	private Vect Fcf;

	private double penalMax;
	private Vect weights;
	public int numContacts;
	public double[] penFactor;
	//public double[] fn_ratio;
	public double[] master_edge_size;
	public double[] fric_coef;
	public int[] type;
	public Node[][] slaveNodes;
	public int[] slaveReg;
	public int[] masterReg;

	public MasterEntity[][]  master_entities;
	boolean stick[], landed_stick[];

	public Element[][] masterElems;

	int[][] normalIndex;
	int[][] u_index;
	int[] numContactingNodes;
	int totalnumContactingNodes;

	boolean[] contacting = null;
	boolean[] just_released = null;
	boolean[] remv = null;


	Vect[][] normals;
	Vect[][] tangentials;

	Model model;
	double pf = 1e8;
	double pft = 1e8;
	double margin = 0.01;

	int n_modifNR0 = 5;
	int aug_itmax0 = 3;
	int nr_itmax0 = 10;
	int nLoads0 = 1;
	double nr_tol0=1e-3;
	double aug_tol0=1e-6;
	
	double nr_tol;
	double aug_tol;
	
	int n_modifNR;
	int aug_itmax;
	int nr_itmax;
	int nLoads;
	
	Vect disp, rhs;
	boolean applyNodal = true;
	boolean initialized = false;
	boolean twice_check =false;
	double fp = 1;
	double fr = 1;
	boolean aug_normal = true;
	boolean aug_tang = true;

	double extention_fact = 0.01;
	double clrFact = 1e-10;
	double aug_disp_tol = 1e-12;
	double gap_tol = 1e-8;
	double reduct = .0;


	double adh = 0e6;
	double adhf = 0e3;

	double relax = .5;;
	
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

		
		if(step==model.nTsteps-1) plot_radial=true;
		
		solver.terminate(false);
		this.model = model;

		direct_slv = new MatSolver();



		fp = 1;
		fr = .01;

		aug_normal = false;
		aug_tang = true;

		applyNodal = true;

		int nmu = 1;
		Vect mus = new Vect(nmu);
		Vect mm = new Vect(nmu);



		KK = null;

		boolean direct = true;


		for (int im = 0; im < nmu; im++) {
		///	disp.zero();
		//	model.setU(disp);
			mus.el[im] =1e3 + .1 * (im);
			for (int contId = 0; contId < numContacts; contId++)
				fric_coef[contId] = mus.el[im];
			


			if(!initialized){
				
				disp = new Vect(model.Ks.nRow);
				if (disp.length > 10000)
					direct = false;
				
			initialize();
			initialized=true;
			
			calcPenaltyFactor();
			}

			rhs =rhs_hat.deepCopy();

			System.out.println(" Contact analysis....");

		//	int[] mm;
			int nout = 0;
			int[] xr_nids1=new int[slaveNodes[0].length];

			for (int k = 0; k <slaveNodes[0].length; k++) {
				Node snode=slaveNodes[0][k];
				if(model.dim==2 ||snode.getCoord(2)<1e-6){
					xr_nids1[nout]=snode.id;
					nout++;
				}
			}

			Vect xr = new Vect(nout);

	
			if (plot_radial)
				for (int k = 0; k <xr.length; k++) {
					int sn=xr_nids1[k];
					Node snode=model.node[sn];
					if(model.dim==2 ||snode.getCoord(2)<1e-6){
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

				for (int cont_iter = 0; cont_iter < aug_itmax; cont_iter++) {

					util.pr("cont_iter: " + cont_iter);

					Vect uaug = disp.deepCopy();

					double er = 1;
					double disp_err = 1;

					boolean zigzag = false;

					for (int nr_iter = 0; nr_iter < nr_itmax; nr_iter++) {

				
						assembleContactMatrices();
						

						addMatrices();
					

						boolean modif_done = false;

						if (/* disp_err<nr_tol && */ totalNRIter > 0 && n_modifNR > 0) {
							if (direct && n_modifNR > 0) {
								KK = Ks.matForm();
								KK.lu();
							}

							modifiedNR(solver, rhs, direct, nr_tol);

							modif_done = true;
						}

						if (modif_done) {
							assembleContactMatrices();
							addMatrices();

						}

						dF = this.calcResidual(Ks, disp, rhs);
		
						er = dF.norm() / rhs.norm();
						nr_err.el[totalNRIter] = er;
						nr_it[totalNRIter] = totalNRIter;

						Vect du = solveLinear(solver, Ks.deepCopy(), dF);

						disp = disp.add(du);

						model.setU(disp);

						if (disp.norm() > 0)
							disp_err = du.norm() / disp.norm();

						util.pr("nr_iter: " + nr_iter + "           nr_err: " + er + "       disp_err: " + disp_err);

						for (int k = 2; k < nr_iter - 1; k += 2) {
							double f = Math.abs(er - nr_err.el[totalNRIter - k]);
							// util.pr("ffff "+f);
							if (f < 1e-2 * nr_tol) {
								zigzag = true;
								break;
							}
						}

						if (zigzag) {
							util.pr("\n ********** NR iteration ended with Zigzag condition! *********\n");
							break;
						}

						totalNRIter++;
						if (er < nr_tol) {
							break;
						}

/*						Vect g = Gc.mul(disp);
						for (int k = 0; k < g.length; k++) {
							if(g.el[k]>0){
								remv[k]=true;
								contacting[k]=false;
							}
							else if(g.el[k]<0) remv[k]=false;
						}*/
						
					}

					double dip_err = disp.sub(uaug).norm() / disp.norm();

					aug_disp_err.el[cont_iter] = dip_err;

					if (dip_err < aug_disp_tol)
						break;

				gap_err.el[cont_iter] = gap.abs().max();

					Vect pgap = gap.deepCopy();
					Vect pslide = slide.deepCopy();// .sub(slide_prev);

					for (int k = 0; k < gap.length; k++) {
						pgap.el[k] *= weights.el[k] * pf;
						pslide.el[k] *= weights.el[k] * pft;
					}

					double ff = 1;
					/// if(gap.abs().max()<gap_tol) ff=0;

					for (int k = 0; k < pgap.length; k++) {

						if (pgap.el[k] >= 0) {
							lamN.el[k] = 0;
							continue;
						}

						double aug = lamN.el[k] + ff * pgap.el[k];


						lamN.el[k] = aug * relax;

					}

					lamT = lamT.add(pslide);

					checkStickSlip();

					if (!aug_normal)
						lamN.zero();

					if (!aug_tang)
						lamT.zero();

					if(Gct!=null)
					aug_N = Gct.mul(lamN);
					
					if(Gcft!=null)
					aug_T = Gcft.mul(lamT);
					

				slide_err.el[cont_iter] = slide.abs().max();

					// xr.show();
					if (plot_radial)
						for (int k = 0; k < xr.length; k++) {
							// int n=2551+k;
							int n = xr_nids[k];

							int index = model.U_unknownIndex[n] - 1;
							if (index < 0)
								continue;

							int com_index = u_index[n][0];
							// urs[cont_iter].el[k]=-Math.abs(u.el[com_index])*1e6;
							if (com_index == -1)
								continue;
							urs[cont_iter].el[k] = disp.el[com_index] * 1e3;
							if (xr.el[k] < 0)
								urs[cont_iter].el[k] *= -1;

							//util.pr(urs[cont_iter].el[k]);
						}

					// if(err.el[cont_iter]<1e-6 && (errf.el[cont_iter]<1e-6))
					// break;

				}
			}


			
			if(step==model.nTsteps-1){
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

			util.pr("slide[m] vs aug_iter");
			slide_err.show("%5.4e");
			// util.plot(errf);

			// ur.show();

			if (plot_radial) {
				double[][] data = new double[xr.length][aug_itmax + 1];
				for (int k = 0; k < xr.length; k++) {

					data[k][0] = xr.el[k];

					for (int i = 0; i < aug_itmax; i++)
						data[k][i + 1] = urs[i].el[k];

				}
				// util.show(data);
				util.plotBunch(data);
				// util.plot(x,ur);
			}
			

			mm.el[im] = disp.abs().max();

		}
		if (nmu > 1)
			util.plot(mus, mm);
		mm.show("%5.4e");

		// s6.show("%5.4e");
		// util.plot(s6);

		// model.setU(new Vect(u.length));

			for (int contId = 0; contId < numContacts; contId++)
				for (int k = 0; k < master_entities[contId].length; k++) {

				
					
					int[] nids = master_entities[contId][k].nodeIds;
		

					Vect normal = normals[contId][k];
					
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
	
		model.writeNodalField(model.resultFolder + "\\master_normal.txt", 2);

		for (int i = 1; i <= model.numberOfNodes; i++) {
			Vect F = model.node[i].Fms;
			if (F != null)
				model.node[i].Fms = null;
		}

		for (int contId = 0; contId < numContacts; contId++)
			for (int k = 0; k < slaveNodes[contId].length; k++) {

				Node node = slaveNodes[contId][k];

				int ind = normalIndex[contId][k];

				if (ind < 0)
					continue;

				Vect normal = normals[contId][ind];// new
													// Vect(-v12.el[1],v12.el[0]).normalized();
				if (normal == null)
					continue;

				model.node[node.id].Fms = normal.times(-1);
				for(int p=0;p<model.dim;p++){
					model.node[node.id].Fms.el[p]+=1e-6*Math.random();
				}
				// node2.u=normal.deepCopy();
			}
		model.writeNodalField(model.resultFolder + "\\slave_normal.txt", 2);

		Fcf=Kcf.smul(disp);
	///	new SpVect(Fc).shownzA();
		if(aug_N.norm()==0)
			model.setU(Fc.times(0).add(Fcf.times(1)).times(-1));
		else
		 model.setU(aug_N.times(0).add(aug_T.times(1)).times(-1));
		for (int i = 1; i <= model.numberOfNodes; i++) {
			// model.node[i].u.el[1]=0;
		}
		
		model.writeNodalField(model.resultFolder + "\\contact_force.txt", -1);
		}

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
			assembleConstraintMats3D();
		///	Gc.shownzA();
		}

		util.pr("totalnumContactingNodes: " + totalnumContactingNodes);

		if (totalnumContactingNodes != 0) {
///Gc.shownzA();
			Gct = Gc.transpose(100);

			Kadh = new SpMat(dof, dof);

			Kc = new SpMat(dof, dof); // Gct*Gc
			
	
			for (int i = 0; i < Gct.nRow; i++) {
				if (Gct.row[i].nzLength > 0) {
					SpVect spv =null;// new SpVect(dof, 100);
					SpVect spvx =null; //new SpVect(dof, 100);

					boolean spv1_filled=false;
					SpVect spv1 =null;// Gct.row[i].deepCopy();

					SpVect spv2 =null;// Gct.row[i].deepCopy();
				

					int kx = 0;

					for (int j = 0; j <= i; j++) {
						if (Gct.row[j].nzLength > 0) {
							
							if(!spv1_filled){
								 spv = new SpVect(dof, 100);
								 spvx = new SpVect(dof, 100);
								 spv1 = Gct.row[i].deepCopy();

								spv2 =Gct.row[i].deepCopy();
								for (int k = 0; k < spv1.nzLength; k++) {
									int ind = spv1.index[k];
								//	util.pr(ind+"  / "+weights.length);
									spv1.el[k] *= weights.el[ind];
								}

								spv1_filled=true;
							}
							double dot = spv1.dot(Gct.row[j]);
							double dotx = spv2.dot(Gct.row[j]);

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
					Kc.row[i] = spv.deepCopy();
					

					spvx.trim(kx);
					Kadh.row[i] = spvx.times(adh);
					}

				}
			}
			
		

			Kc.times(pf);


			Gcft = Gcf.transpose(100);
			G_stkt = G_stk.transpose(100);

			Kcf = new SpMat(dof, dof); // Gct*Gc

			// === adhisve tang
			Kadhf = new SpMat(dof, dof);
			SpMatAsym Gcfadht = Gcfadh.transpose(100);

			for (int i = 0; i < Gcfadht.nRow; i++) {
				if (Gcfadht.row[i].nzLength > 0) {

					SpVect spv = new SpVect(dof, 100);

					SpVect spv1 = Gcfadht.row[i].deepCopy();

					int kx = 0;
					for (int j = 0; j <= i; j++) {
						if (Gcfadht.row[j].nzLength > 0) {

							double dot = spv1.dot(Gcfadht.row[j]);

							if (dot == 0)
								continue;

							spv.index[kx] = j;
							spv.el[kx++] = dot;

						}
					}

					spv.trim(kx);
					Kadhf.row[i] = spv.times(adhf);

				}
			}

			// ====

			for (int i = 0; i < G_stkt.nRow; i++) {
				if (G_stkt.row[i].nzLength > 0) {
					SpVect spv = new SpVect(dof, 100);
					SpVect spv1 = G_stkt.row[i].deepCopy();

					for (int k = 0; k < spv1.nzLength; k++) {
						int ind = spv1.index[k];
						spv1.el[k] *= weights.el[ind];
					}

					int kx = 0;
					for (int j = 0; j <= i; j++) {
						if (G_stkt.row[j].nzLength > 0) {

							double dot = spv1.dot(G_stkt.row[j]);

							if (dot == 0)
								continue;
							spv.index[kx] = j;
							spv.el[kx++] = dot;
						}
					}
					spv.trim(kx);
					Kcf.row[i] = spv.deepCopy();

				}
			}

			Kcf.times(pft);

		}
	}

	private void addMatrices() {

		if (Kc != null) {

			Ks = K_hat.addGeneral(Kc);
			
			Ks = Ks.addGeneral(Kadh.addGeneral(Kadhf));
		} else {
			Ks = K_hat.deepCopy();
		}
		
	}

	private void modifiedNR(SpMatSolver solver, Vect bU1, boolean direct, double tol) {

		double er1 = 1;
		Vect uc[] =new Vect[n_modifNR];
		Vect errs=new Vect().ones(n_modifNR).times(-1);
		for (int sb = 0; sb < n_modifNR; sb++) {

			/*
			 * if(sb%5==1){ assembleContactMatrices(); addMatrices();
			 * 
			 * }else
			 */
			if(model.dim==2)
			updateGap(true);
			else
			{
				obtain_node_node3D();
				assembleConstraintMats3D();
		
			}
	

			Vect dF1 = calcResidual(Ks, disp, bU1);

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

	}

	private void obtain_node_node() {

		// slide_prev.zero();
		gap.zero();

		// weights.zero();
		// weightsf.zero();

		numContactingNodes = new int[numContacts];

		totalnumContactingNodes = 0;

		int nnSize = model.numberOfNodes + 1;
		for (int contId = 0; contId < numContacts; contId++) {
			

			
			double mu = fric_coef[contId];
			
			for (int i = 0; i < slaveNodes[contId].length; i++) {
				Node node = slaveNodes[contId][i];
				int sn = node.id;
				
				if(remv[sn]){
					remv[sn]=false;
					util.pr("------------------- 888 ");
					continue;
				}

				node_node.row[sn] = new SpVect(nnSize);

				Vect u = node.u;
				Vect v = node.getCoord().add(u);

				if (master_edge_size[contId] == 0) {
					double length = 0;
					for (int k = 0; k < master_entities[contId].length; k++) {

						int[] nids=master_entities[contId][k].nodeIds;
						
						Node node1 = model.node[nids[0]];
						Node node2 = model.node[nids[1]];

						Vect v12 = node1.getCoord().add(node2.getCoord());

						length += v12.norm();
					}
					master_edge_size[contId] = length;
				}

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids=master_entities[contId][k].nodeIds;
					
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

						contacting[sn] = false;

						continue;
					}

					Element elem = masterElems[contId][k];
					int[] vn = elem.getVertNumb();
					for (int j = 0; j < vn.length; j++) {
						if (vn[j] != mn1 && vn[j] != mn2) {
							Node node3 = model.node[vn[j]];

							Vect v3 = node3.getCoord().add(node3.u);
							Vect v13 = v3.sub(v1);

							Vect cross1 = edgeDir.v3().cross(v13.v3());
							Vect cross2 = edgeDir.v3().cross(cross1.v3());
							normals[contId][k] = cross2.normalized().v2();
							break;
						}
					}

					Vect normal = normals[contId][k];

					double pen = v1v.dot(normal);

					if (pen < -100 * edgeLength) {
						// weakenning.el[p]=0;

						contacting[sn] = false;

						continue;
					}

					// if(!gradualSeperation){
					if (pen > clrFact * master_edge_size[contId]) {



					//	if (!twice_check || !just_released[sn]){
				
					//	if(contacting[sn]) just_released[sn]=true;

						contacting[sn] = false;

						continue;
					//	}
						
					}


					normalIndex[contId][i] = k;

					double beta = v1v.dot(edgeDir) / edgeLength;
					double alpha = 1 - beta;

					node_node.row[sn] = new SpVect(nnSize, 2);
					node_node.row[sn].index[0] = node1.id;
					node_node.row[sn].index[1] = node2.id;
					node_node.row[sn].el[0] = alpha;
					node_node.row[sn].el[1] = beta;

	

					Vect deltaDisp = u.sub(u1.times(alpha).add(u2.times(beta)));
				//	Vect deltaDisp = v.sub(v1.times(alpha).add(v2.times(beta)));
					
					Vect vr=new Vect(u.length);
					Vect v1r=new Vect(u.length);
					Vect v2r=new Vect(u.length);
					if(u_index[sn][0]>=0)
					vr.el[0]=	ref_stick.el[u_index[sn][0]];
					if(u_index[sn][1]>=0)
						vr.el[1]=	ref_stick.el[u_index[sn][1]];
					
					if(u_index[mn1][0]>=0)
					v1r.el[0]=	ref_stick.el[u_index[mn1][0]];
					if(u_index[mn1][1]>=0)
						v1r.el[1]=	ref_stick.el[u_index[mn1][1]];
					
					if(u_index[mn2][0]>=0)
						v2r.el[0]=	ref_stick.el[u_index[mn2][0]];
						if(u_index[mn2][1]>=0)
							v2r.el[1]=	ref_stick.el[u_index[mn2][1]];
						
				
					Vect deltaDispRef = vr.sub(v1r.times(alpha).add(v2r.times(beta)));
					
				//	deltaDisp=deltaDisp.sub(deltaDispRef);
					double proj = deltaDisp.dot(normal);

					Vect disp_tang = deltaDisp.sub(normal.times(proj));
					
				//	disp_tang.show("%10.5e");

					// util.pr("-----< "+disp_tang.norm());

					Vect tang = disp_tang.normalized();
					// if(tang.norm()==0)
					tang = edgeDir.deepCopy();

					tangentials[contId][k] = tang.deepCopy();
	
					gap.el[sn] = pen;
					slide.el[sn] = deltaDisp.dot(tang);



					contacting[sn] = true;

					totalnumContactingNodes++;

					numContactingNodes[contId]++;
					break;//
				}

			}
		}

	//	resetFreedNodes();

		countStickSlip();
		/// new SpVect(weights).shownzA();
		for (int contId = 0; contId < numContacts; contId++)
			util.pr((contId + 1) + ", num slave nodes: " + slaveNodes[contId].length + ",  numContactingNodes: "
					+ numContactingNodes[contId]);

		// new SpVect(weights).shownzA();
	}

	private void obtain_node_node3D() {

		// slide_prev.zero();
		gap.zero();

		numContactingNodes = new int[numContacts];

		totalnumContactingNodes = 0;

		int nnSize = model.numberOfNodes + 1;
		for (int contId = 0; contId < numContacts; contId++) {
			double mu = fric_coef[contId];
			for (int i = 0; i < slaveNodes[contId].length; i++) {
				Node node = slaveNodes[contId][i];
				int sn = node.id;

				node_node.row[sn] = new SpVect(nnSize);

				Vect u = node.u;
				Vect v = node.getCoord().add(u);

				if (master_edge_size[contId] == 0) {
					double length = 0;
					for (int k = 0; k < master_entities[contId].length; k++) {

						int[] nids=master_entities[contId][k].nodeIds;
						
						Node node1 = model.node[nids[0]];
						Node node2 = model.node[nids[1]];
						Node node3 = model.node[nids[2]];
						Vect v12 = node1.getCoord().add(node3.getCoord());

						double len=v12.norm();
						
						master_entities[contId][k].length=len;
						length += len;
						
					}
					master_edge_size[contId] = length;
				}

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids=master_entities[contId][k].nodeIds;

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
				
					normal=new Vect(0,1,0);
					normals[contId][k] = normal.deepCopy();

				//	 normal.hshow();

					Vect v1v=v.sub(v1);
					Vect cv=v.sub(center);
					double pen = cv.dot(normal);


					if (pen > clrFact * master_edge_size[contId]) {

						contacting[sn] = false;

						continue;
					}

					normalIndex[contId][i] = k;

					
					
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
					
					 double[] ww1=new double[4];
					 double uu1=.5;
					 double vv1=.5;
					 ww1[0]=(1-uu1)*(1-vv1);
					 ww1[1]=(uu1)*(1-vv1);
					 ww1[2]=(uu1*vv1);
					 ww1[3]=(1-uu1)*(vv1);
					
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


					node_node.row[sn] = new SpVect(nnSize, 4);
					for (int m = 0; m < 4; m++) {
						node_node.row[sn].index[m] = nids[m];
						///if (m == nnn)
						node_node.row[sn].el[m] = ww[m];

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

					gap.el[sn] = pen;
					

					Vect deltaDisp = u.deepCopy();
					for(int m=0;m<nnc;m++){
						deltaDisp=deltaDisp.sub(um[m].times(ww[m]));
						
					}
					//u.sub(u1.times(alpha).add(u2.times(beta)));
					/// deltaDisp.times(1e9).hshow();
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

					tangentials[contId][k] = tang.deepCopy();

					//double sld=(disp_tang.dot(tang));
					double sld=Math.abs(disp_tang.dot(tang));

					slide.el[sn]=sld;
					
					contacting[sn] = true;

					totalnumContactingNodes++;

					numContactingNodes[contId]++;
					break;//
				}

			}
		}

	 //node_node.shownzA();

		//resetFreedNodes();

		countStickSlip();
		/// new SpVect(weights).shownzA();
		for (int contId = 0; contId < numContacts; contId++)
			util.pr((contId + 1) + ", num slave nodes: " + slaveNodes[contId].length + ",  numContactingNodes: "
					+ numContactingNodes[contId]);

		//
	///new SpVect(gap).shownzA();
	}

	private void updateGap(boolean allowSep) {

		gap.zero();

		slide.zero();

		for (int contId = 0; contId < numContacts; contId++) {

			for (int i = 0; i < slaveNodes[contId].length; i++) {
				Node node = slaveNodes[contId][i];
				int sn = node.id;

				if (!contacting[sn])
					continue;

				Vect u = node.u;
				Vect v = node.getCoord().add(u);

				if (node_node.row[sn].index == null)
					continue;

				int mn1 = node_node.row[sn].index[0];
				int mn2 = node_node.row[sn].index[1];

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

				int nrmIndex = normalIndex[contId][i];

				Vect normal = null;
				Element elem = masterElems[contId][nrmIndex];
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

				if (allowSep && pen > clrFact * master_edge_size[contId]) {

					continue;
				}

				if (pen < -100 * edgeLength) {

					continue;
				}

				double alpha = node_node.row[sn].el[0];
				double beta = node_node.row[sn].el[1];

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

				gap.el[sn] = pen;

				slide.el[sn] = deltaDisp.dot(tang);

			}
		}

		resetFreedNodes();

	}

	private void assembleConstraintMats() {

		int dof = model.Ks.nRow;
		int nRows=model.numberOfNodes + 1;
		Gc = new SpMatAsym(nRows, dof);
		Gcf = new SpMatAsym(nRows, dof);
		Gcfadh = new SpMatAsym(nRows, dof);

		G_stk = new SpMatAsym(nRows, dof);

		for (int contId = 0; contId < numContacts; contId++)
			for (int i = 0; i < slaveNodes[contId].length; i++) {

				Node node = slaveNodes[contId][i];
				int sn = node.id;

				if (node_node.row[sn].nzLength > 0) {

					int index = model.U_unknownIndex[sn] - 1;
					if (index < 0)
						continue;

					int px = u_index[sn][0];
					int py = u_index[sn][1];
					int pz = -1;
					if (model.dim == 3)
						pz = u_index[sn][2];

	

					Gc.row[sn] = new SpVect(dof);

					Vect normal = normals[contId][normalIndex[contId][i]];

				
					int[] nids=master_entities[contId][normalIndex[contId][i]].nodeIds;

					int mn1 = nids[0];
					int mn2 = nids[1];

					int p1x = u_index[mn1][0];
					int p1y = u_index[mn1][1];

					int p2x = u_index[mn2][0];
					int p2y = u_index[mn2][1];

					double alpha = node_node.row[sn].el[0];
					double beta = node_node.row[sn].el[1];

					Gc.row[sn] = new SpVect(dof, 6);

					int kx = 0;
					if (px != -1) {
						Gc.row[sn].index[kx] = px;
						Gc.row[sn].el[kx++] = normal.el[0];
					}
					if (py != -1) {
						Gc.row[sn].index[kx] = py;
						Gc.row[sn].el[kx++] = normal.el[1];
					}

					// util.pr(weights.el[com_index]);

					if (p1x >= 0) {
						Gc.row[sn].index[kx] = p1x;
						Gc.row[sn].el[kx++] = -alpha * normal.el[0];
					}
					if (p1y >= 0) {
						Gc.row[sn].index[kx] = p1y;
						;
						Gc.row[sn].el[kx++] = -alpha * normal.el[1];
					}

					if (p2x >= 0) {
						Gc.row[sn].index[kx] = p2x;
						Gc.row[sn].el[kx++] = -beta * normal.el[0];
					}
					if (p2y >= 0) {
						Gc.row[sn].index[kx] = p2y;
						;
						Gc.row[sn].el[kx++] = -beta * normal.el[1];
					}
					Gc.row[sn].sortAndTrim(kx);
					

					// Vect tang=new Vect(-normal.el[1],normal.el[0]);

					Vect tang = tangentials[contId][normalIndex[contId][i]].deepCopy();

					// addhesive tangential{
					// ====
					Gcfadh.row[sn] = new SpVect(dof, 6);
					kx = 0;
					if (px != -1) {
						Gcfadh.row[sn].index[kx] = px;
						Gcfadh.row[sn].el[kx++] = tang.el[0];
					}
					if (p1y != -1) {
						Gcfadh.row[sn].index[kx] = py;
						Gcfadh.row[sn].el[kx++] = tang.el[1];
					}

					if (p1x >= 0) {
						Gcfadh.row[sn].index[kx] = p1x;
						Gcfadh.row[sn].el[kx++] = -alpha * tang.el[0];
					}
					if (p1y >= 0) {
						Gcfadh.row[sn].index[kx] = p1y;
						;
						Gcfadh.row[sn].el[kx++] = -alpha * tang.el[1];
					}

					if (p2x >= 0) {
						Gcfadh.row[sn].index[kx] = p2x;
						Gcfadh.row[sn].el[kx++] = -beta * tang.el[0];
					}
					if (p2y >= 0) {
						Gcfadh.row[sn].index[kx] = p2y;
						;
						Gcfadh.row[sn].el[kx++] = -beta * tang.el[1];
					}

					// ===

					if (fric_coef[contId] == 0) {
						tang.zero();
					}
					Gcf.row[sn] = new SpVect(dof, 6);
					kx = 0;
					if (px != -1) {
						Gcf.row[sn].index[kx] = px;
						Gcf.row[sn].el[kx++] = tang.el[0];
					}
					if (p1y != -1) {
						Gcf.row[sn].index[kx] = py;
						Gcf.row[sn].el[kx++] = tang.el[1];
					}

					if (p1x >= 0) {
						Gcf.row[sn].index[kx] = p1x;
						Gcf.row[sn].el[kx++] = -alpha * tang.el[0];
					}
					if (p1y >= 0) {
						Gcf.row[sn].index[kx] = p1y;
						;
						Gcf.row[sn].el[kx++] = -alpha * tang.el[1];
					}

					if (p2x >= 0) {
						Gcf.row[sn].index[kx] = p2x;
						Gcf.row[sn].el[kx++] = -beta * tang.el[0];
					}
					if (p2y >= 0) {
						Gcf.row[sn].index[kx] = p2y;
						;
						Gcf.row[sn].el[kx++] = -beta * tang.el[1];
					}

					if (stick[sn]) {
						G_stk.row[sn] = Gcf.row[sn].deepCopy();

					} else {
						G_stk.row[sn] = Gcf.row[sn].times(reduct);
					}

				}
			}

		// G_stk.shownzA();

	}

	private void assembleConstraintMats3D() {

		int dof = model.Ks.nRow;

		Gc = new SpMatAsym(model.numberOfNodes + 1, dof);
		Gcf = new SpMatAsym(model.numberOfNodes + 1, dof);
		Gcfadh = new SpMatAsym(model.numberOfNodes + 1, dof);

		G_stk = new SpMatAsym(model.numberOfNodes + 1, dof);

		for (int contId = 0; contId < numContacts; contId++)
			for (int i = 0; i < slaveNodes[contId].length; i++) {

				Node node = slaveNodes[contId][i];
				int sn = node.id;

				if (node_node.row[sn].nzLength > 0) {

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


					Vect normal = normals[contId][normalIndex[contId][i]];

				
					int[] nids = master_entities[contId][normalIndex[contId][i]].nodeIds;
							

					int p1x = u_index[nids[0]][0];
					int p1y = u_index[nids[0]][1];
					int p1z = u_index[nids[0]][2];

					int p2x = u_index[nids[1]][0];
					int p2y = u_index[nids[1]][1];
					int p2z = u_index[nids[1]][2];

					int p3x = u_index[nids[2]][0];
					int p3y = u_index[nids[2]][1];
					int p3z = u_index[nids[2]][2];

					int p4x = u_index[nids[3]][0];
					int p4y = u_index[nids[3]][1];
					int p4z = u_index[nids[3]][2];

					double alpha = node_node.row[sn].el[0];
					double beta = node_node.row[sn].el[1];
					double gamma = node_node.row[sn].el[2];
					double zeta = node_node.row[sn].el[3];

					Gc.row[sn] = new SpVect(dof, 15);

					int kx = 0;
					if (px != -1) {
						Gc.row[sn].index[kx] = px;
						Gc.row[sn].el[kx++] = normal.el[0];
					}
					if (py != -1) {
						Gc.row[sn].index[kx] = py;
						Gc.row[sn].el[kx++] = normal.el[1];
					}
					if (pz != -1) {
						Gc.row[sn].index[kx] = pz;
						Gc.row[sn].el[kx++] = normal.el[2];
					}
				
					// util.pr(weights.el[com_index]);

					if (p1x >= 0) {
						Gc.row[sn].index[kx] = p1x;
						Gc.row[sn].el[kx++] = -alpha * normal.el[0];
					}
					if (p1y >= 0) {
						Gc.row[sn].index[kx] = p1y;
						;
						Gc.row[sn].el[kx++] = -alpha * normal.el[1];

					}
					if (p1z >= 0) {
						Gc.row[sn].index[kx] = p1z;
						;
						Gc.row[sn].el[kx++] = -alpha * normal.el[2];

					}

					if (p2x >= 0) {
						Gc.row[sn].index[kx] = p2x;
						Gc.row[sn].el[kx++] = -beta * normal.el[0];
					}
					if (p2y >= 0) {
						Gc.row[sn].index[kx] = p2y;
						;
						Gc.row[sn].el[kx++] = -beta * normal.el[1];
					}
					if (p2z >= 0) {
						Gc.row[sn].index[kx] = p2z;
					
						Gc.row[sn].el[kx++] = -beta * normal.el[2];
					}

					if (p3x >= 0) {
						Gc.row[sn].index[kx] = p3x;
						Gc.row[sn].el[kx++] = -gamma * normal.el[0];
					}
					if (p3y >= 0) {
						Gc.row[sn].index[kx] = p3y;
						;
						Gc.row[sn].el[kx++] = -gamma * normal.el[1];
					}
					if (p3z >= 0) {
						Gc.row[sn].index[kx] = p3z;
						;
						Gc.row[sn].el[kx++] = -gamma * normal.el[2];
					}

					if (p4x >= 0) {
						Gc.row[sn].index[kx] = p4x;
						Gc.row[sn].el[kx++] = -zeta * normal.el[0];
					}
					if (p4y >= 0) {
						Gc.row[sn].index[kx] = p4y;
						;
						Gc.row[sn].el[kx++] = -zeta * normal.el[1];
					}
					if (p4z >= 0) {
						Gc.row[sn].index[kx] = p4z;
						;
						Gc.row[sn].el[kx++] = -zeta * normal.el[2];
					}
				
					Gc.row[sn].sortAndTrim(kx);

					// Vect tang=new Vect(-normal.el[1],normal.el[0]);

					Vect tang = tangentials[contId][normalIndex[contId][i]].deepCopy();

				

					if (fric_coef[contId] == 0) {
						tang.zero();
					}
					Gcf.row[sn] = new SpVect(dof, 15);

					 kx = 0;
					if (px != -1) {
						Gcf.row[sn].index[kx] = px;
						Gcf.row[sn].el[kx++] = tang.el[0];
					}
					if (py != -1) {
						Gcf.row[sn].index[kx] = py;
						Gcf.row[sn].el[kx++] = tang.el[1];
					}
					if (pz != -1) {
						Gcf.row[sn].index[kx] = pz;
						Gcf.row[sn].el[kx++] = tang.el[2];
					}
				
					

					if (p1x >= 0) {
						Gcf.row[sn].index[kx] = p1x;
						Gcf.row[sn].el[kx++] = -alpha * tang.el[0];
					}
					if (p1y >= 0) {
						Gcf.row[sn].index[kx] = p1y;
						;
						Gcf.row[sn].el[kx++] = -alpha * tang.el[1];

					}
					if (p1z >= 0) {
						Gcf.row[sn].index[kx] = p1z;
						;
						Gcf.row[sn].el[kx++] = -alpha * tang.el[2];

					}

					if (p2x >= 0) {
						Gcf.row[sn].index[kx] = p2x;
						Gcf.row[sn].el[kx++] = -beta * tang.el[0];
					}
					if (p2y >= 0) {
						Gcf.row[sn].index[kx] = p2y;
						;
						Gcf.row[sn].el[kx++] = -beta * tang.el[1];
					}
					if (p2z >= 0) {
						Gcf.row[sn].index[kx] = p2z;
					
						Gcf.row[sn].el[kx++] = -beta * tang.el[2];
					}

					if (p3x >= 0) {
						Gcf.row[sn].index[kx] = p3x;
						Gcf.row[sn].el[kx++] = -gamma * tang.el[0];
					}
					if (p3y >= 0) {
						Gcf.row[sn].index[kx] = p3y;
						;
						Gcf.row[sn].el[kx++] = -gamma * tang.el[1];
					}
					if (p3z >= 0) {
						Gcf.row[sn].index[kx] = p3z;
						;
						Gcf.row[sn].el[kx++] = -gamma * tang.el[2];
					}

					if (p4x >= 0) {
						Gcf.row[sn].index[kx] = p4x;
						Gcf.row[sn].el[kx++] = -zeta * tang.el[0];
					}
					if (p4y >= 0) {
						Gcf.row[sn].index[kx] = p4y;
						;
						Gcf.row[sn].el[kx++] = -zeta * tang.el[1];
					}
					if (p4z >= 0) {
						Gcf.row[sn].index[kx] = p4z;
						;
						Gcf.row[sn].el[kx++] = -zeta * tang.el[2];
					}
					
					if (stick[sn]) {
						G_stk.row[sn] = Gcf.row[sn].deepCopy();

					} else {
						G_stk.row[sn] = Gcf.row[sn].times(reduct);
					}
				}
			}

		// G_stk.shownzA();
	}

	private void checkStickSlip() {

		for (int contId = 0; contId < numContacts; contId++) {
			double mu = this.fric_coef[contId];
			if (mu == 0)
				continue;
			for (int i = 0; i < slaveNodes[contId].length; i++) {

				Node node = slaveNodes[contId][i];
				int sn = node.id;

				if (!contacting[sn]) {
					stick[sn] = false;
					landed_stick[sn] = false;
					continue;
				}

				int index = model.U_unknownIndex[sn] - 1;
				if (index < 0)
					continue;

				if (node_node.row[sn].nzLength > 0) {
					if (lamN.el[sn] < 0) {
						double abs_lamT = Math.abs(lamT.el[sn]);
						double muFn = mu * Math.abs(lamN.el[sn]);
						// util.ph("node "+sn+ " p "+p+ " lamT: ");
						// util.pr(lamT.el[p]," %4.4e");
						// util.ph("node "+sn+ " lamN: ");
						// util.pr(lamN.el[p]," %4.4e");
						if (abs_lamT > muFn * (1 + margin)) {
							if (lamT.el[sn] > 0)
								lamT.el[sn] = muFn;
							else
								lamT.el[sn] = -muFn;
							stick[sn] = false;
							landed_stick[sn] = false;
							slide.el[sn] = 0;
						} else {
							if (!stick[sn]) {
								stick[sn] = true;

								landed_stick[sn] = true;
								slide.el[sn] = 0;
							}
						}
					} else {
						stick[sn] = false;
						landed_stick[sn] = false;
						lamT.el[sn] = 0;
						slide.el[sn] = 0;
						contacting[sn] = false;
					}
				} else {
					stick[sn] = false;
					landed_stick[sn] = false;
					lamT.el[sn] = 0;
					slide.el[sn] = 0;
					contacting[sn] = false;
				}
			}

		}
		
		for (int contId = 0; contId < numContacts; contId++) {
			double mu = this.fric_coef[contId];
			if (mu == 0)
				continue;
			for (int i = 0; i < slaveNodes[contId].length; i++) {

				Node node = slaveNodes[contId][i];
				int sn = node.id;


			if (landed_stick[sn]) {
				

				int px = u_index[sn][0];
				int py = u_index[sn][1];
				int pz = -1;
				if (model.dim == 3)
					pz = u_index[sn][2];

				int p = px;
				if (p == -1)
					p = py;
				if (p == -1)
					p = pz;


				int[] nids=master_entities[contId][normalIndex[contId][i]].nodeIds;

				int mn1 = nids[0];
				int mn2 = nids[1];

				int p1x = u_index[mn1][0];
				int p1y = u_index[mn1][1];

				int p2x = u_index[mn2][0];
				int p2y = u_index[mn2][1];
				
				Vect u = model.node[sn].getU();
				if (px >= 0)
					ref_stick.el[px] = u.el[0];
				if (py >= 0)
					ref_stick.el[py] = u.el[1];

				u = model.node[mn1].getU();
				if (p1x >= 0)
					ref_stick.el[p1x] = u.el[0];
				if (p1y >= 0)
					ref_stick.el[p1y] = u.el[1];

				u = model.node[mn2].getU();
				if (p2x >= 0)
					ref_stick.el[p2x] = u.el[0];
				if (p2y >= 0)
					ref_stick.el[p2y] = u.el[1];

				landed_stick[sn] = false;

			}
		}
		
			}
		countStickSlip();

		/*
		 * int px=0,qx=0;
		 * 
		 * for(int i=0;i<Gc.nRow;i++){ if(Gc.row[i].nzLength>0){ px++; } }
		 * for(int i=0;i<G_stk.nRow;i++){ if(G_stk.row[i].nzLength>0){ qx++; } }
		 * 
		 * util.pr(", px: "+px+",  qx: "+qx);
		 */
		/// G_stkt.shownzA();

		// G_stk.shownzA();

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

	private Vect calcResidual(SpMat Ks, Vect u, Vect b) {

		Vect Fint = K_hat.smul(u);

		Vect dF = b.sub(Fint);


		Vect pg = gap.times(pf);
	//	new SpVect(pg).shownzA();
		for (int k = 0; k < gap.length; k++) {
			pg.el[k] *= weights.el[k];
		}

		if(Gct!=null)
		Fc = Gct.mul(pg);


		Vect ps = slide.times(pft);
		for (int k = 0; k < gap.length; k++) {
			ps.el[k] *= weights.el[k];

		}
	//	G_stk.shownzA();

		if(G_stkt!=null)
		Fcf = G_stkt.mul(ps);

		dF = dF.sub(Fc).sub(Fcf); //

		dF = dF.sub(aug_N).sub(aug_T);

		return dF;

	}

	public void readContact(Loader loader, BufferedReader br, Model model) throws IOException {

		model.setEdge();
		

		double minEdgeLenghth = model.minEdgeLength;
		util.pr(" minEdgeLength ------------------------- " + minEdgeLenghth);
		double clearFact = 1e-4;
		
		String line = loader.getNextDataLine(br," / * NUM. CONTACT * NR_ITR. * AUG_ITER  * NLOADS  * MNR /");

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


		numContacts = numCont;
		
		nr_tol=nr_tol0;
		aug_tol=aug_tol0;
		
		if(integs.length>1){
			
			 line = loader.getNextDataLine(br," / * NR_TOL  *  AUG_TOL * /");
			 double[] tols=loader.getTabedData(line);
			 nr_tol=tols[0];
			 if(integs.length>2)
			 aug_tol=tols[1];
		}

		slaveNodes = new Node[numCont][];
		
		slaveReg=new int[numCont];
		masterReg=new int[numCont];
		
		master_entities = new MasterEntity[numCont][];

		penFactor = new double[numCont];
		fric_coef = new double[numCont];

		master_edge_size = new double[numCont];

		for (int i = 0; i < numCont; i++) {
			line = loader.getNextDataLine(br," /* PENALTY FACTOR */");
			penFactor[i] = loader.getScalarData(line);
			line = loader.getNextDataLine(br," /* COEF OF FRICTION */");
			fric_coef[i] = loader.getScalarData(line);
			line =loader.getNextDataLine(br);

			int type = 0;
			if (!line.contains("slav"))
				type = loader.getIntData(line);

			if (type == 0) {
				line = loader.getNextDataLine(br," /* NUM SLAVE NODES */");

				int ns = loader.getIntData(line);

				slaveNodes[i] = new Node[ns];

				for (int k = 0; k < ns; k++) {
					line = br.readLine();

					int sn = loader.getIntData(line);
					slaveNodes[i][k] = model.node[sn];
				}
			} else {
				if (model.dim == 3) {
					//util.pr("This format of setting contact not ready yet.");
					line = loader.getNextDataLine(br," /* REGION ID */");
					int nreg = loader.getIntData(line);
					
					line = loader.getNextDataLine(br," /* CORENR 1 */");

					int n1 = loader.getIntData(line);

					line = loader.getNextDataLine(br," /* CORENR 2 */");

					int n2 = loader.getIntData(line);
					
					line = loader.getNextDataLine(br," /* CORENR 3 */");

					int n3 = loader.getIntData(line);
					
					line = loader.getNextDataLine(br,"/* CORENR 4 */");

					int n4 = loader.getIntData(line);

					Node node1 = model.node[n1];
					Node node2 = model.node[n2];
					Node node3 = model.node[n3];
					Node node4 = model.node[n4];
					
			
			

					Vect v1 = node1.getCoord();
					Vect v2 = node2.getCoord();
					Vect v3 = node3.getCoord();
					Vect v4 = node4.getCoord();
					
		
					Vect v13=v3.sub(v1);
					Vect v24=v4.sub(v2);
					
					v1 = v1.add(v13.times(-extention_fact));
					v3 = v3.add(v13.times(extention_fact));
					v2 = v2.add(v24.times(-extention_fact));
					v4 = v4.add(v24.times(extention_fact));
					
					Vect v12=v2.sub(v1);
					Vect v23=v3.sub(v2);
					Vect v34=v4.sub(v3);
					Vect v41=v1.sub(v4);

						
					int ns = 0;
					int[] nnr = model.getRegNodes(nreg);
					int[] temp = new int[nnr.length];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						Vect v = model.node[n].getCoord();
						Vect v1v = v.sub(v1);
						Vect v2v = v.sub(v2);
						Vect v3v = v.sub(v3);
						Vect v4v = v.sub(v4);
						
				
						
						Vect cross[]=new Vect[4];
						
						cross[0]=v12.cross(v1v);
						cross[1]=v23.cross(v2v);
						cross[2]=v34.cross(v3v);
						cross[3]=v41.cross(v4v);
						
						
						boolean infront=true;
						for(int j=0;j<4;j++){
							for(int m=j+1;m<4;m++){
								double dot=cross[j].dot(cross[m]);
								if(dot<0) {
									infront=false;
									break;
								}
								}
							if(!infront) break;
							}
						
						if(!infront) continue;
						
						Vect normal=v13.cross(v23);
						double dist=Math.abs(v1v.dot(normal));
						
						if(dist>clearFact * minEdgeLenghth) continue;


								temp[ns] = n;
								ns++;
							
						

					}
					
					util.pr("contact "+i+":  num. slave nodes ------------------------- " + ns);
					
					slaveNodes[i] = new Node[ns];

					for (int k = 0; k < ns; k++) {

						int sn = temp[k];
						slaveNodes[i][k] = model.node[sn];
						
						//util.pr(sn);
					}
				//	util.pr(ns);


				} else {
					line = br.readLine();
					int nreg = loader.getIntData(line);
					
					slaveReg[i] = nreg;
				
					
					line = br.readLine();

					int n1 = loader.getIntData(line);

					line = br.readLine();

					int n2 = loader.getIntData(line);

					Node node1 = model.node[n1];
					Node node2 = model.node[n2];

					Vect v1 = node1.getCoord();
					Vect v2 = node2.getCoord();
					Vect edgeDir = v2.sub(v1).normalized();
					int ns = 0;

					int[] nnr = model.getRegNodes(nreg);
					int[] temp = new int[nnr.length];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						Vect v = model.node[n].getCoord();
						Vect vv1 = v.sub(v1);
						Vect vv2 = v.sub(v2);
						double dot = vv1.dot(vv2);

						if (dot <= 0) {

							double proj = vv1.dot(edgeDir);
							double vv1n = vv1.norm();
							if (Math.abs(proj - vv1n) < 1e-4 * minEdgeLenghth) {
								// util.pr(dot);

								temp[ns] = n;
								ns++;
							}
						}

					}
					slaveNodes[i] = new Node[ns];

					for (int k = 0; k < ns; k++) {

						int sn = temp[k];
						slaveNodes[i][k] = model.node[sn];
					}

				}
			}

			line = br.readLine();
			type = 0;
			if (!line.contains("mast"))
				type = loader.getIntData(line);
			if (type == 0) {
				line = br.readLine();

				int nm = loader.getIntData(line);

		
				master_entities[i] = new MasterEntity[nm];


				for (int k = 0; k < nm; k++) {
					line = br.readLine();

					int[] nn = loader.getCSInt(line);
				//	if (model.dim == 2) {
			
						Node node1=model.node[nn[0]];
						Node node2=model.node[nn[1]];
						master_entities[i][k]=new MasterEntity(2);
						master_entities[i][k].nodeIds[0]=nn[0];
						master_entities[i][k].nodeIds[1]=nn[1];
						master_entities[i][k].length= node1.getCoord().sub(node2.getCoord()).norm();
	/*				} else {

						
						if (nn.length == 3) {
							masterFacets[i][k] = new Element("triangle");
							masterFacets[i][k].setVertNumb(nn);

						} else if (nn.length == 4) {
							masterFacets[i][k] = new Element("quad");
							masterFacets[i][k].setVertNumb(nn);

						}
					}*/
				}
			} else {
				if (model.dim == 3) {
				//	util.pr("This format of setting contact not ready yet.");
					
					byte[][] arr0 = { { 0, 1,2,3 }, { 4,7,6,5 }, { 0,4,5,1 },{2,6,7,3},{0,3,7,4 }, {1,5,6,2}};
				//	byte[][] arr1 = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };

					byte[][] edgeLocalNodes = null;
					//if (model.elCode == 0) {
						edgeLocalNodes = arr0;
					//	;
				//	} else if (model.elCode == 1) {
					//	edgeLocalNodes = arr1;
				//	}
						
		
						line = loader.getNextDataLine(br," /* REGION ID */");
						int nreg = loader.getIntData(line);
						
						masterReg[i] = nreg;
						
						line = loader.getNextDataLine(br," /* CORENR 1 */");

						int n1 = loader.getIntData(line);

						line = loader.getNextDataLine(br," /* CORENR 2 */");

						int n2 = loader.getIntData(line);
						
						line = loader.getNextDataLine(br," /* CORENR 3 */");

						int n3 = loader.getIntData(line);
						
						line = loader.getNextDataLine(br,"/* CORENR 4 */");

						int n4 = loader.getIntData(line);
						Node node1 = model.node[n1];
						Node node2 = model.node[n2];
						Node node3 = model.node[n3];
						Node node4 = model.node[n4];
						
				

						Vect v1 = node1.getCoord();
						Vect v2 = node2.getCoord();
						Vect v3 = node3.getCoord();
						Vect v4 = node4.getCoord();
						
						Vect v13=v3.sub(v1);
						Vect v24=v4.sub(v2);
						
						v1 = v1.add(v13.times(-extention_fact));
						v3 = v3.add(v13.times(extention_fact));
						v2 = v2.add(v24.times(-extention_fact));
						v4 = v4.add(v24.times(extention_fact));
						
						Vect v12=v2.sub(v1);
						Vect v23=v3.sub(v2);
						Vect v34=v4.sub(v3);
						Vect v41=v1.sub(v4);

						
						double totalArea=v12.cross(v23).norm()/2+v34.cross(v41).norm()/2;

						int[] nnr = model.getRegNodes(nreg);
						boolean[] nc = new boolean[model.numberOfNodes + 1];
						for (int k = 0; k < nnr.length; k++) {
							int n = nnr[k];
							nc[n] = true;
						}
						
						int nm=0;
						
						int[][] tempEd = new int[model.numberOfEdges * edgeLocalNodes.length][4];
						for (int ie = model.region[nreg].getFirstEl(); ie <= model.region[nreg].getLastEl(); ie++) {
							int[] vertNumb = model.element[ie].getVertNumb();
							
							
							
				

							for (int j = 0; j < edgeLocalNodes.length; j++) {
								
								Vect cent=new Vect(3);
								for(int k=0;k<4;k++){
									int nx = vertNumb[edgeLocalNodes[j][k]];
									Vect v = model.node[nx].getCoord();
									cent=cent.add(v.times(.25));
								}

								Vect v = cent;
								Vect v1v = v.sub(v1);
								Vect v2v = v.sub(v2);
								Vect v3v = v.sub(v3);
								Vect v4v = v.sub(v4);
									
									Vect cross[]=new Vect[4];
									
									cross[0]=v12.cross(v1v);
									cross[1]=v23.cross(v2v);
									cross[2]=v34.cross(v3v);
									cross[3]=v41.cross(v4v);
									
									
									boolean infront=true;
									for(int k=0;k<4;k++){
										for(int m=k+1;m<4;m++){
											double dot=cross[k].dot(cross[m]);
											if(dot<0) {
												infront=false;
												break;
											}
											}
										if(!infront) break;
										}
									
									if(!infront) continue;
									
									Vect normal=v13.cross(v23);
									double dist=Math.abs(v1v.dot(normal));
									
									if(dist>1e-4  * minEdgeLenghth) continue;
					
								for(int p=0;p<4;p++){
									int n12 = vertNumb[edgeLocalNodes[j][p]];
						
										tempEd[nm][p] = n12;
									
									}
							//	util.hshow(tempEd[nm]);
								nm++;
							}

						}

						util.pr("contact "+i+":  num. master facets ------------------------- " + nm);
						master_entities[i] = new MasterEntity[nm];

						for (int k = 0; k < nm; k++) {

							master_entities[i][k] = new MasterEntity(4);
							for(int p=0;p<4;p++){
								master_entities[i][k].nodeIds[p]=tempEd[k][p];
							}
						}
							
						
						
				} else {
					line = br.readLine();
					int nreg = loader.getIntData(line);
					
					masterReg[i] = nreg;
					
					line = br.readLine();
					int n1 = loader.getIntData(line);
					line = br.readLine();
					int n2 = loader.getIntData(line);
					Node node1 = model.node[n1];
					Node node2 = model.node[n2];
					Vect v1 = node1.getCoord();
					Vect v2 = node2.getCoord();
					Vect edgeDir = v2.sub(v1).normalized();
					int nm = 0;

					int[] nnr = model.getRegNodes(nreg);
					boolean[] nc = new boolean[model.numberOfNodes + 1];
					for (int k = 0; k < nnr.length; k++) {
						int n = nnr[k];
						nc[n] = true;
					}

					byte[][] arr0 = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
					byte[][] arr1 = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };

					byte[][] edgeLocalNodes = null;
					if (model.elCode == 0) {
						edgeLocalNodes = arr0;
						
					} else if (model.elCode == 1) {
						edgeLocalNodes = arr1;
					}

					int[][] tempEd = new int[model.numberOfEdges * edgeLocalNodes.length][2];
					for (int ie = model.region[nreg].getFirstEl(); ie <= model.region[nreg].getLastEl(); ie++) {
						int[] vertNumb = model.element[ie].getVertNumb();

						for (int j = 0; j < edgeLocalNodes.length; j++) {
							int n11 = vertNumb[edgeLocalNodes[j][0]];
							int n12 = vertNumb[edgeLocalNodes[j][1]];

							Vect v = model.node[n11].getCoord().add(model.node[n12].getCoord()).times(0.5);
							Vect vv1 = v.sub(v1);
							Vect vv2 = v.sub(v2);
							double dot = vv1.dot(vv2);

							if (dot <= 0) {

								double proj = vv1.dot(edgeDir);
								double vv1n = vv1.norm();
								if (Math.abs(proj - vv1n) < clearFact * minEdgeLenghth) {
									// util.pr(dot);

									tempEd[nm][0] = n11;
									tempEd[nm][1] = n12;
									nm++;
								}
							}
						}

					}

					master_entities[i] = new MasterEntity[nm];

					for (int k = 0; k < nm; k++) {

						Node node11 = model.node[tempEd[k][0]];
						Node node12 = model.node[tempEd[k][1]];

						master_entities[i][k]=new MasterEntity(2);
						master_entities[i][k].nodeIds[0]=tempEd[k][0];
						master_entities[i][k].nodeIds[1]=tempEd[k][1];
				
						master_entities[i][k].length= node11.getCoord().sub(node12.getCoord()).norm();

					}
				}

			}

		}
	
	boolean node_duplic=false;
	if(node_duplic){
		
		for (int contId = 0; contId < numContacts; contId++) {
			Node[] sns = slaveNodes[contId];
			MasterEntity[]  med = master_entities[contId];
			
			boolean[] coincid=new boolean[model.numberOfNodes+1];

			for (int i = 0; i < coincid.length; i++) 
				coincid[i]=false;
			
			for (int i = 0; i < sns.length; i++) {
				int sn = slaveNodes[contId][i].id;
				for (int j = 0; j < med.length; j++) {
					int n1=med[j].nodeIds[0];
					int n2=med[j].nodeIds[1];
				///	util.pr(sn+" "+n1+"  "+n2);
					if(sn==n1){
						coincid[sn]=true;
						break;
					}
					else if(sn==n2){
						coincid[sn]=true;	
						break;
					}
				
				}
		}
			
			int[] map=new int[model.numberOfNodes+1];
			for (int i = 0; i < map.length; i++)
				map[i]=0;
			
			int extera_nodes=0;
			for (int i = 0; i < coincid.length; i++) {
				if(coincid[i]) {
					extera_nodes++;
					map[i]=model.numberOfNodes+extera_nodes;
				//	util.pr(i+" ---- "+map[i]);
				}
			}

			Vect[] dup_coord=new Vect[extera_nodes];
			int ix=0;

			for (int i = 0; i < map.length; i++) {
				if(map[i]>0) dup_coord[ix++]=model.node[i].getCoord();
			}
			
			int nRegions=model.numberOfRegions;
			int nElements=model.numberOfElements;
			int nNodes=model.numberOfNodes+extera_nodes;
			Model md1=new Model(nRegions,nElements,nNodes,model.elType);
	
			int n1=1,N;
			for(int i=1;i<=md1.numberOfRegions;i++){

				md1.region[i].setFirstEl(model.region[i].getFirstEl());
				md1.region[i].setLastEl(model.region[i].getLastEl());
				md1.region[i].setName(model.region[i].getName());
				md1.region[i].setMaterial(model.region[i].getMaterial());


			}

			for(int i=1;i<=model.numberOfNodes;i++)	
				md1.node[i].setCoord(model.node[i].getCoord());
			
			for(int i=1;i<=extera_nodes;i++)	{
				md1.node[i+model.numberOfNodes].setCoord(dup_coord[i-1]);
			}

			for(int ir=1;ir<=md1.numberOfRegions;ir++)		{

		
				boolean on_master=masterReg[contId]==ir;
				
				for(int i=md1.region[ir].getFirstEl();i<=md1.region[ir].getLastEl();i++)	{
					int[] vn=model.element[i].getVertNumb();
						
					for(int k=0;k<md1.nElVert;k++){
						if(map[vn[k]]==0 )
							md1.element[i].setVertNumb(k,vn[k]);
						else{
							if(on_master)
								md1.element[i].setVertNumb(k,vn[k]);
							else
								md1.element[i].setVertNumb(k,map[vn[k]]);
						
						}
					}

				
			}
			}



			md1.scaleFactor=model.scaleFactor;
			
			String folder=new File(model.meshFilePath).getParentFile().getPath();
			md1.meshFilePath=model.meshFilePath;
			
			model=md1.deepCopy();
			model.meshFilePath=md1.meshFilePath;
		for (int contId2 = 0; contId2 < numContacts; contId2++) {
			for (int i = 0; i <slaveNodes[contId2].length; i++) {
				int sn = slaveNodes[contId2][i].id;
			
				if(map[sn]>0){ slaveNodes[contId2][i]=model.node[map[sn]];
				//util.pr(sn+"  "+slaveNodes[contId][i].id);
				}
			}
			if(contId2!=contId)
			for (int i = 0; i <master_entities[contId2].length; i++) {
				
				int mn1 = master_entities[contId2][i].nodeIds[0];
				int mn2 = master_entities[contId2][i].nodeIds[1];
			
				if(map[mn1]>0){
					master_entities[contId2][i].nodeIds[0]=map[mn1];
				}
				if(map[mn2]>0){
					master_entities[contId2][i].nodeIds[1]=map[mn2];
				}
			}
		}
			String file = folder + "//duplicated.txt";
			if(contId==numCont-1){
				int[] nns=model.getRegNodes(7);
				for (int i = 0; i < nns.length; i++) {
					int n = nns[i];
					Vect v=model.node[n].getCoord();
					v.el[0]+=-.001;
					v.el[1]+=.001;
					
					model.node[n].setCoord(v);
				
					
				}
		/*		nns=model.getRegNodes(2);
				for (int i = 0; i < nns.length; i++) {
					int n = nns[i];
					Vect v=model.node[n].getCoord();
					v.el[1]+=-.05;
					
					model.node[n].setCoord(v);
				
					
				}*/
			}
			model.writeMesh(file);
		

	}
		
		
		
	}

	}

	private void resetFreedNodes() {
		for (int contId = 0; contId < numContacts; contId++) {
			for (int i = 0; i < slaveNodes[contId].length; i++) {
				int sn = slaveNodes[contId][i].id;
				if (!contacting[sn]) {
					stick[sn] = false;
					landed_stick[sn] = false;

					lamN.el[sn] = 0;
					lamT.el[sn] = 0;

				}
			}
		}

	}

	private void calcPenaltyFactor() {

		penalMax = 0;

		for (int contId = 0; contId < numContacts; contId++)
			for (int i = 0; i < slaveNodes[contId].length; i++) {
				int sn = slaveNodes[contId][i].id;

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
					val1 = fp * model.Ks.row[xind].el[model.Ks.row[xind].nzLength - 1];

				if (yind != -1)
					val2 = fp * model.Ks.row[yind].el[model.Ks.row[yind].nzLength - 1];
				if (zind != -1)
					val3 = fp * model.Ks.row[zind].el[model.Ks.row[zind].nzLength - 1];

				if (val1 > penalMax)
					penalMax = val1;
				if (val2 > penalMax)
					penalMax = val2;
				if (val3 > penalMax)
					penalMax = val3;

				// double max=(val1+val2)/2;
			//	double max = Math.max(Math.max(val1, val2), val3);
			double max=Math.sqrt(val1*val1+val2*val2+val3*val3);
				// util.pr(max);
				// max=1e10;
				//
	

				if (applyNodal)
					weights.el[sn] = penFactor[contId] * max;
				else
					weights.el[sn] = penFactor[contId];

				// weightsf.el[p]=weights.el[p]*fn_ratio[contId];

				////// weights.el[p]*=1./slaveNodes[contId].length;

			}

		// weights0=weights.deepCopy();

		pf = penalMax;
		pft = fr * pf;

		util.pr("pf :" + pf);
		util.pr("pft :" + pft);
		if (applyNodal) {
			pf = 1e0;
			pft = fr * pf;
		}
	}

	private void countStickSlip() {

		for (int contId = 0; contId < slaveNodes.length; contId++) {
			int nstk = 0;
			int nslip = 0;
			for (int i = 0; i < slaveNodes[contId].length; i++) {
				Node node = slaveNodes[contId][i];
				int sn = node.id;

				if (contacting[sn]) {
	/*				int px = u_index[sn][0];
					int py = u_index[sn][1];

					int p = px;
					if (p == -1)
						p = py;*/
					// double abs_lamT=Math.abs(lamT.el[p]);
					// double muFn=mu*Math.abs(lamN.el[p]);
					// util.pr("lamT.el[p] "+lamT.el[p] +" lamN.el[p]
					// "+lamN.el[p]);

					if (stick[sn])
						nstk++;
					else
						nslip++;

					// util.pr("node "+sn+" stick "+stick[sn]);
					/// util.pr("u "+model.node[sn].getU(p)+" uref
					// "+ref_stick.el[p]);
				}
				// else
				// util.pr("node "+sn+" free ");
			}
			util.pr((contId + 1) + ", stick: " + nstk + ",  slip: " + nslip);
		}
	}

	private void initialize() {

		int[][] ne = new int[model.numberOfNodes + 1][20];
		int[] nz = new int[model.numberOfNodes + 1];

		for (int i = 1; i <= model.numberOfElements; i++) {
			int[] vn = model.element[i].getVertNumb();
			for (int j = 0; j < vn.length; j++) {
				ne[vn[j]][nz[vn[j]]++] = i;

			}
		}

		masterElems = new Element[numContacts][];
		for (int contId = 0; contId < numContacts; contId++) {
			if (model.dim == 2) {
				masterElems[contId] = new Element[master_entities[contId].length];

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids=master_entities[contId][k].nodeIds;
					
					Node node1 = model.node[nids[0]];
					Node node2 = model.node[nids[1]];

					int ie = 0;
					for (int j = 0; j < nz[node1.id]; j++) {
						for (int p = 0; p < nz[node2.id]; p++) {

							if (ne[node1.id][j] == ne[node2.id][p]) {
								ie = ne[node1.id][j];
								break;
							}
							if (ie > 0)
								break;
						}
					}

					if (ie > 0)
						masterElems[contId][k] = model.element[ie];
					else
						util.pr("master edge ( " + node1.id + ", " + node2.id + " ) belongs to no element.");

				}
			} else {

				masterElems[contId] = new Element[master_entities[contId].length];

				for (int k = 0; k < master_entities[contId].length; k++) {

					int[] nids = master_entities[contId][k].nodeIds;

					int ie = 0;
					for (int j = 0; j < nz[nids[0]]; j++) {
						for (int p = 0; p < nz[nids[1]]; p++) {
							for (int q = 0; q < nz[nids[2]]; q++) {

								if (ne[nids[0]][j] == ne[nids[1]][p] && ne[nids[0]][j] == ne[nids[2]][q]) {
									ie = ne[nids[0]][j];

									break;

								}
								if (ie > 0)
									break;
							}
							if (ie > 0)
								break;
						}
					}

					if (ie > 0) {

						masterElems[contId][k] = model.element[ie];
					} else
						util.pr("master facet ( " + nids[0] + ", " + nids[1] + ", " + nids[2] + ", " + nids[3]
								+ " ) belongs to no element.");

				}
			}

			// for(int k=0;k<masterEdges[contId].length;k++)
			// util.hshow(masterElems[contId][k].getVertNumb());
		}

		landed_stick = new boolean[model.numberOfNodes + 1];

		stick = new boolean[model.numberOfNodes + 1];
		contacting = new boolean[model.numberOfNodes + 1];
		just_released = new boolean[model.numberOfNodes + 1];
		remv = new boolean[model.numberOfNodes + 1];
		
		for (int i = 1; i <= model.numberOfNodes; i++) {
			stick[i] = true;
			landed_stick[i] = true;
		}

		u_index = new int[model.numberOfNodes + 1][model.dim];
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
						ix++;
					}
				}

			}
		}


		ref_stick = new Vect(model.Ks.nRow);

		lamN = new Vect(model.numberOfNodes + 1);
		lamT = new Vect(model.numberOfNodes + 1);

		aug_N = new Vect(model.Ks.nRow);
		aug_T = new Vect(model.Ks.nRow);

		weights = new Vect(model.numberOfNodes + 1);// .ones(model.Ks.nRow);
		// weightsf=new Vect(model.Ks.nRow);
		// weakenning= new Vect(model.Ks.nRow).ones(model.Ks.nRow);
		// weakenningf= new Vect(model.Ks.nRow).ones(model.Ks.nRow);

		gap = new Vect(model.numberOfNodes + 1);
		slide = new Vect(model.numberOfNodes + 1);

		// slide_prev=new Vect(model.Ks.nRow);

		Fc = new Vect(model.Ks.nRow);
		Fcf = new Vect(model.Ks.nRow);

		normalIndex = new int[numContacts][];

		normals = new Vect[numContacts][];
		tangentials = new Vect[numContacts][];
		
		type=new int[numContacts];

		for (int contId = 0; contId < numContacts; contId++) {
			
			type[contId]=0;
			
			int numSn = slaveNodes[contId].length;
			int numMed = 0;
	
			numMed = master_entities[contId].length;
			normalIndex[contId] = new int[numSn];
			normals[contId] = new Vect[numMed];
			tangentials[contId] = new Vect[numMed];

		}

		int nnSize = model.numberOfNodes + 1;
		node_node = new SpMatAsym(nnSize, nnSize);
		
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

		double loadFactor0 = 1000;

		if (model.dim == 3)
			loadFactor0 = 1;
		
		if (model.centrigForce) {
			model.setNodalMass();
			double rpm = 7000;
			double rps = rpm / 30 * PI;
			double omeag2 = rps * rps;
			for (int i = 1; i <= model.numberOfNodes; i++) {
				Vect v = model.node[i].getCoord();
				double m = model.node[i].getNodalMass();
				Vect F = v.times(omeag2 * m/loadFactor0);
				model.node[i].setF(F);
			}

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

	
	public static void main(String[] args) {

		new Main();
	}

}

