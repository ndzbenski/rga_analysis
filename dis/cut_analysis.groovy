import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// usage:
// rungroovy xB_analysis.groovy <filename.hipo> <energy>

double en = Double.parseDouble(args[1]);
double enmax = en+0.1; //GeV
double thetamax = 40;  //degrees
double phimax = 180;   //degrees
double vzmax = 50;
double wmin = 0;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

HipoDataSource reader = new HipoDataSource();

// The reconstructed 1D histos 
H1F theta_hist = new H1F("theta", "theta", 50, 0, thetamax+5);
theta_hist.setTitleX("theta [deg]");

H1F phi_hist = new H1F("phi", "phi", 100, -phimax, phimax);
phi_hist.setTitleX("phi [deg]");

H1F momentum = new H1F("momentum", "momentum", 20, 0, 11);
momentum.setTitleX("momentum [GeV]");

H1F W_hist = new H1F("W", "W", 100, wmin, wmax+0.5);
W_hist.setTitleX("W [GeV]");

H1F Q2_hist = new H1F("Q2", "Q2", 50, 0, 13);
Q2_hist.setTitleX("Q^2 [GeV^2]");

H1F Eprime_hist = new H1F("Eprime", "E'", 50, 0, 13);
Eprime_hist.setTitleX("E' [GeV]");

H1F xB_hist = new H1F("xB", "xB", 20, 0, 1);
xB_hist.setTitleX("xB");

// Create multiple histograms for each Q2 bin 
HashMap<Integer,H1F> xB_histmap = new HashMap<Integer,H1F>();
for(int i = 0; i < 8; i++){xB_histmap.put(i,new H1F("xB", 100, 0, 0.9));}

H1F xsection_hist = new H1F("xsection", "Generating Cross Section (#sigma)", 100, 0, 1);
xsection_hist.setTitleX("#sigma_{gen}");


// The 1D histos cuts
H1F theta_hist_cut = new H1F("theta_cut", "theta_cut", 50, 0, thetamax+5);
theta_hist_cut.setTitleX("theta [deg]");

H1F phi_hist_cut = new H1F("phi_cut", "phi_cut", 100, -phimax, phimax);
phi_hist_cut.setTitleX("phi [deg]");

H1F mom_hist_cut = new H1F("momentum_cut", "momentum_cut", 20, 0, 11);
mom_hist_cut.setTitleX("p [GeV]");

H1F W_hist_cut = new H1F("W_cut", "W_cut", 100, wmin, wmax+0.5);
W_hist_cut.setTitleX("W_cut [GeV]");

H1F Q2_hist_cut = new H1F("Q2_cut", "Q2_cut", 100, 0, 13);
Q2_hist_cut.setTitleX("Q^2_cut [GeV^2]");

H1F xB_hist_cut = new H1F("xB_cut", "xB_cut", 20, 0, 1);
xB_hist_cut.setTitleX("xB_cut");


// 2D Histos
H2F Q2_vs_W = new H2F("Q2_vs_W", "Q2 vs W", 100, wmin, wmax+0.5, 100, 0.0, 13);
Q2_vs_W.setTitleX("W [GeV]");
Q2_vs_W.setTitleY("Q^2 [GeV^2]");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E' vs Theta", 100, 5, thetamax+5, 100, 0, enmax);
E_vs_Theta.setTitleX("Theta [deg]");
E_vs_Theta.setTitleY("E' [GeV]");

H2F Q2_vs_xB = new H2F("Q2_vs_xB", "Q2 vs xB", 100, 0, 1, 100, 0, 13);
Q2_vs_xB.setTitleX("xB");
Q2_vs_xB.setTitleY("Q^2 [GeV^2]");

H2F W_vs_xB = new H2F("W_vs_xB", "W vs xB", 100, 0, 0.81, 100, wmin, wmax+0.5);
W_vs_xB.setTitleX("xB");
W_vs_xB.setTitleY("W [GeV]");

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 150, wmin, wmax, 500, -phimax, phimax);
Phi_vs_W.setTitleX("W [GeV]");
Phi_vs_W.setTitleY("Phi [deg]");

H2F xsect_vs_xB = new H2F("xsect_vs_xB", "generating #sigma vs xB", 150, 0.0, 0.81, 150, -1, 150);
xsect_vs_xB.setTitleX("xB");
xsect_vs_xB.setTitleY("#sigma");


// For pre-fiducial cuts
H2F Cal_y_vs_x_precut = new H2F("Cal_y_vs_x_precut", "Cal_y_vs_x_precut", 100, -450,450, 100, -450,450);
Cal_y_vs_x_precut.setTitleX("X [cm]");
Cal_y_vs_x_precut.setTitleY("Y [cm]");

H1F Cal_lu_precut = new H1F("Cal_lu", "Cal_lu_precut", 50, 0, 450);
Cal_lu_precut.setTitleX("Cal_lu_precut [cm]");
H1F Cal_lv_precut = new H1F("Cal_lv", "Cal_lv_precut", 50, 0, 450);
Cal_lv_precut.setTitleX("Cal_lv_precut [cm]");
H1F Cal_lw_precut = new H1F("Cal_lw", "Cal_lw_precut", 50, 0, 450);
Cal_lw_precut.setTitleX("Cal_lw_precut [cm]");

// fiducial cuts
H1F Cal_lu = new H1F("Cal_lu", "Cal_lu", 50, 0, 450);
Cal_lu.setTitleX("Cal_lu_cut [cm]");
H1F Cal_lv = new H1F("Cal_lv", "Cal_lv", 50, 0, 450);
Cal_lv.setTitleX("Cal_lv_cut [cm]");
H1F Cal_lw = new H1F("Cal_lw", "Cal_lw", 50, 0, 450);
Cal_lw.setTitleX("Cal_lw_cut [cm]");

H2F Cal_y_vs_x = new H2F("Cal_y_vs_x", "Cal_y_vs_x", 100, -450,450, 100, -450, 450);
Cal_y_vs_x.setTitleX("X [cm]");
Cal_y_vs_x.setTitleY("Y [cm]");


double e_mass = 0.000511;
double p_mass = 0.93827203;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);


// create for loop for all files
// open cooked.files
// read in line by line
// for each line, open and run analysis
// close file
new File('/work/clas12/nated/dis.cooked/', args[0]).eachLine { line ->
    reader.open(line);
    
    double emax = 0;
    phimax = 0;
    thetamax = 0;
    vzmax = 0;
    int counter = 0;
    byte sector = 0;
    int cal_row = 0;
    int dc_row = 0;
    float Q2_bin_min = 0;
    float Q2_bin_max = 0;
    
    // generated data variables
    float weight = 0;
    
    // reconstructed data variables
    int pid = 0;
    byte q = 0;
    float px = 0;
    float py = 0;
    float pz = 0;
    float beta =0;
    float mom = 0;
    double phi = 0;
    double theta =  0;
    float vz = 0;
    double Q2 = 0; 
    double W =0;
    double E_prime = 0;
    double xB = 0;
                
                
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        double tot_xsect = 0;
        
        // get reconstructed data
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj") && event.hasBank("MC::Event")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            
            //System.out.println("number of bank_rec rows: " + bank_rec.rows() + ", # of gen bank rows: " + bank_gen.rows() );
                
            for (int k = 0; k < bank_rec.rows(); k++) {
    
                pid = bank_rec.getInt("pid", k);
                q = bank_rec.getByte("charge", k);
                px = bank_rec.getFloat("px", k);
                py = bank_rec.getFloat("py", k);
                pz = bank_rec.getFloat("pz", k);
                beta = bank_rec.getFloat("beta", k);
    
                mom = (float) Math.sqrt(px * px + py * py + pz * pz);
                phi = Math.atan2((double) py,(double) px);
                theta = Math.acos((double) pz/(double) mom);
                
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                vz = bank_rec.getFloat("vz", k);
    
                // pick electrons
                //if (pid != 11) continue;
                if(q != -1) continue;
                
                Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
                LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
                e_vec_prime.setVectM(e_vec_3, e_mass);
    
                // ---------------- Cut on E' and theta -------------------
                //if(e_vec_prime.e() < 0.1 * en){continue;}  //cut below 10% beam
                //if(theta < 5 || theta > 40){continue;}     //cut outside of 5 and 40 degrees for FD
                // --------------------------------------------------------

                cal_row = cal_cut_row(event, k);
                //System.out.println(j + " " + bank_cal.rows());
                if(cal_row != -1){
                    sector = bank_cal.getByte("sector",cal_row);
                    
                    float x_cal = bank_cal.getFloat("x",cal_row);
                    float y_cal = bank_cal.getFloat("y",cal_row);
                    
                    float lu = bank_cal.getFloat("lu",cal_row);
                    float lv = bank_cal.getFloat("lv",cal_row);
                    float lw = bank_cal.getFloat("lw",cal_row);
                    
                    Cal_y_vs_x_precut.fill(x_cal,y_cal);
                    Cal_lu_precut.fill(lu);
                    Cal_lv_precut.fill(lv);
                    Cal_lw_precut.fill(lw);
                    
                    if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
                        Cal_lu.fill(lu);
                        Cal_lv.fill(lv);
                        Cal_lw.fill(lw);
                        Cal_y_vs_x.fill(x_cal,y_cal);
                    }
                }
                
                dc_row = dc_cut_row(event, k);
                if(dc_row != -1){
                    float x_dc = bank_traj.getFloat("x",dc_row);
                    float y_dc = bank_traj.getFloat("y",dc_row);
                    float z_dc = bank_traj.getFloat("z",dc_row);
                    
                    double pos = Math.sqrt(x_dc*x_dc + y_dc*y_dc + z_dc*z_dc);
                    double theta_dc = Math.acos((double) z_dc/ pos);
                    double phi_dc = Math.atan2((double) y_dc,(double) x_dc);
                    
                    theta_dc *= 180/Math.PI;
                    phi_dc *= 180/Math.PI;
                    
                }

                momentum.fill(mom);
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                Q2 = -q_vec.mass2(); //-q^2
                
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                W = w_vec.mass();
                //double W = p_mass*p_mass + 2.0*p_mass*(en-E_prime) -Q2
    
                //double E_prime = Q2/(4.0*en*(Math.sin(theta*Math.PI/360.0)));
                E_prime = e_vec_prime.e();
                xB = Q2/(2.0*p_mass*(en-E_prime));
                
                xB_hist.fill(xB);
                
                // Fill histos
                Q2_hist.fill(Q2);
                Eprime_hist.fill(E_prime);
                
                // xB histograms binned in 1 GeV^2 Q2
                for(int j = 0; j < 8; j++){
                    Q2_bin_min = 0;
                    Q2_bin_max = 1*j + 1.0;
                    if(Q2 < Q2_bin_max) {
                        xB_histmap.get(j).fill(xB);
                        
                    }
                    else {continue;}
                }
                
                W_hist.fill(W);
                xsection_hist.fill(weight);
                theta_hist.fill(theta);
                phi_hist.fill(phi);
                
                // Calculate max values of each param
                if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} 
                if(theta > thetamax){thetamax = theta;}
                if(phi > phimax){phimax = phi;}
                if(vz > vzmax){vzmax = vz;}
                
                // begin cuts
                if (theta < 5 || theta > 40) {continue;}  
                if (W < 2) {continue;}
                if (E_prime < 0.1*en) {continue;}
                
                theta_hist_cut.fill(theta);
                phi_hist_cut.fill(phi);
                mom_hist_cut.fill(mom);
                    
                W_hist_cut.fill(W);
                Q2_hist_cut.fill(Q2);
                xB_hist_cut.fill(xB);
                
                Q2_vs_W.fill(W,Q2);
                Phi_vs_W.fill(W,phi);
                E_vs_Theta.fill(theta,e_vec_prime.e()); // check this with calculated E'
                Q2_vs_xB.fill(xB,Q2);
                W_vs_xB.fill(xB,W);
                xsect_vs_xB.fill(xB,weight);
                
            } // end for
        } // end if    
    } // end while
    
} // end open file

boolean dc_cut(float X, float Y, int S){
    boolean result= false;
    if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
    && (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
    
    return result;
}

int cal_cut_row(DataEvent event, int row){
    DataBank bank_cal = event.getBank("RECHB::Calorimeter");
    
    int row_index = 0;
    int cal_row_match = -1;
    
    for(int j = 0; j < bank_cal.rows(); j++){
        row_index = bank_cal.getInt("pindex",j);
        
        if(row_index == row){
            cal_row_match = j;
            break;
        }
    }
    return cal_row_match;
}

int dc_cut_row(DataEvent event, int row){
    DataBank bank_traj = event.getBank("REC::Traj");
    int row_index = 0;
    int det_id = 0;
    int cal_row_match = -1;
    for(int j = 0; j < bank_traj.rows(); j++){
        row_index = bank_traj.getInt("pindex",j);
        det_id = bank_traj.getInt("detId",j);
        
        if(row_index == row && det_id == 6){
            cal_row_match = j;
            break;
        }
    }
    return cal_row_match;
}

TCanvas can_1d_a = new TCanvas("can", 1100, 600);
can_1d_a.setTitle("Theta, phi & mom cuts");
can_1d_a.divide(3,2);
can_1d_a.cd(0);
can_1d_a.draw(theta_hist); 
can_1d_a.cd(1);
can_1d_a.draw(phi_hist);
can_1d_a.cd(2);
can_1d_a.draw(momentum);
can_1d_a.cd(3);
can_1d_a.draw(theta_hist_cut);
can_1d_a.cd(4);
can_1d_a.draw(phi_hist_cut);
can_1d_a.cd(5);
can_1d_a.draw(mom_hist_cut);

//for(int k =0; k < 8; k++){ 
//    xB_histmap.get(k).setLineColor(k);
//    can_1d.draw(xB_histmap.get(k),"same");
    //can.getPad().setLegend(true);
    //can.getPad().setLegendPosition(20, 20);
    //legend.AddEntry(xB_histmap.get(k),"Q2 < " + 1*k+1,"l");
//}
can_1d_a.save("figs/cuts/1D_angle_cuts.png");

TCanvas can_1d_b = new TCanvas("can", 1100, 600);
can_1d_b.setTitle("W, Q2, xB resolutions");
can_1d_b.divide(3,2);
can_1d_b.cd(0);
can_1d_b.draw(W_hist); 
can_1d_b.cd(1);
can_1d_b.draw(Q2_hist);
can_1d_b.cd(2);
can_1d_b.draw(xB_hist);
can_1d_b.cd(3);
can_1d_b.draw(W_hist_cut);
can_1d_b.cd(4);
can_1d_b.draw(Q2_hist_cut);
can_1d_b.cd(5);
can_1d_b.draw(xB_hist_cut);
can_1d_b.save("figs/cuts/1D_kin_cuts.png");

TCanvas can_2d_a = new TCanvas("can", 1100, 600);
can_2d_a.divide(2,2);
can_2d_a.cd(0);
can_2d_a.draw(Q2_vs_W);
can_2d_a.cd(1);
can_2d_a.draw(E_vs_Theta);
can_2d_a.cd(2);
can_2d_a.draw(Q2_vs_xB);
can_2d_a.cd(3);
can_2d_a.draw(W_vs_xB);
can_2d_a.save("figs/cuts/2d_spectra_cuts.png");

TCanvas can_2d_b = new TCanvas("can", 1100, 600);
can_2d_b.divide(4,2);
can_2d_b.cd(0);
can_2d_b.draw(Cal_lu_precut);
can_2d_b.cd(1);
can_2d_b.draw(Cal_lv_precut);
can_2d_b.cd(2);
can_2d_b.draw(Cal_lw_precut);
can_2d_b.cd(3);
can_2d_b.draw(Cal_y_vs_x_precut);
can_2d_b.cd(4);
can_2d_b.draw(Cal_lu);
can_2d_b.cd(5);
can_2d_b.draw(Cal_lv);
can_2d_b.cd(6);
can_2d_b.draw(Cal_lw);
can_2d_b.cd(7);
can_2d_b.draw(Cal_y_vs_x);
can_2d_b.save("figs/cuts/2d_ecal.png");

/*
HashMap<Integer,TCanvas> canvasmap = new HashMap<Integer,TCanvas>();
for(int i : histmap.keySet()){
canvasmap.put(i, new TCanvas("can", 800,600));
canvasmap.get(i).draw(histmap.get(i));
canvasmap.get(i).save("phivstheta" + i + ".png");
}*/