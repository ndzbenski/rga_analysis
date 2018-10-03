import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
//import org.jlab.groot.data.legend;
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

// The generated 1D histos 
H1F theta_hist_gen = new H1F("theta_gen", "theta_gen", 500, 0, thetamax+5);
theta_hist_gen.setTitleX("theta_gen [deg]");

H1F phi_hist_gen = new H1F("phi_gen", "phi_gen", 500, -phimax, phimax);
phi_hist_gen.setTitleX("phi_gen [deg]");

H1F mom_hist_gen = new H1F("gen momentum", "gen momentum", 500, 0, 11);
mom_hist_gen.setTitleX("momentum_gen [GeV]");

H1F W_hist_gen = new H1F("W_gen", "W_gen", 500, wmin, wmax+0.5);
W_hist_gen.setTitleX("W_gen [GeV]");

H1F Q2_hist_gen = new H1F("Q2_gen", "Q2_gen", 50, 0, 13);
Q2_hist_gen.setTitleX("Q^2 gen [GeV^2]");

H1F xB_hist_gen = new H1F("xB_gen", "xB_gen", 100, 0, 1);
xB_hist_gen.setTitleX("xB_gen");


// The reconstructed 1D histos 
H1F theta_hist = new H1F("theta", "theta", 500, 0, thetamax+5);
theta_hist.setTitleX("theta [deg]");

H1F phi_hist = new H1F("phi", "phi", 500, -phimax, phimax);
phi_hist.setTitleX("phi [deg]");

H1F momentum = new H1F("momentum", "momentum", 500, 0, 11);
momentum.setTitleX("momentum [GeV]");

H1F W_hist = new H1F("W", "W", 500, wmin, wmax+0.5);
W_hist.setTitleX("W [GeV]");

H1F Q2_hist = new H1F("Q2", "Q2", 50, 0, 13);
Q2_hist.setTitleX("Q^2 [GeV^2]");

H1F Eprime_hist = new H1F("Eprime", "E'", 50, 0, 13);
Eprime_hist.setTitleX("E' [GeV]");

H1F xB_hist = new H1F("xB", "xB", 100, 0, 1);
xB_hist.setTitleX("xB");

// Create multiple histograms for each Q2 bin 
HashMap<Integer,H1F> xB_histmap = new HashMap<Integer,H1F>();
for(int i = 0; i < 8; i++){xB_histmap.put(i,new H1F("xB", 100, 0, 0.9));}

H1F xsection_hist = new H1F("xsection", "Generating Cross Section (#sigma)", 100, 0, 1);
xsection_hist.setTitleX("#sigma_{gen}");


// The resolution 1D histos 
H1F theta_hist_res = new H1F("theta_res", "theta_res", 500, -2, 2);
theta_hist_res.setTitleX("#Delta theta [deg]");

H1F phi_hist_res = new H1F("phi_res", "phi_res", 500, -10, 10);
phi_hist_res.setTitleX("#Delta phi [deg]");

H1F mom_hist_res = new H1F("momentum_res", "momentum_res", 500, -1, 1);
mom_hist_res.setTitleX("#Delta p [GeV]");

H1F W_hist_res = new H1F("W_res", "W_res", 500, -1, 1);
W_hist_res.setTitleX("W_res [GeV]");

H1F Q2_hist_res = new H1F("Q2_res", "Q2_res", 500, -1, 1);
Q2_hist_res.setTitleX("Q^2_res [GeV^2]");

H1F xB_hist_res = new H1F("xB_res", "xB_res", 100, -1, 1);
xB_hist_res.setTitleX("xB_res");


// 2D Histos
H2F Q2_vs_W = new H2F("Q2_vs_W", "Q2 vs W", 500, wmin, wmax+0.5, 500, 0.0, 13);
Q2_vs_W.setTitleX("W [GeV]");
Q2_vs_W.setTitleY("Q^2 [GeV^2]");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E' vs Theta", 500, 5, thetamax+5, 500, 0, enmax);
E_vs_Theta.setTitleX("Theta [deg]");
E_vs_Theta.setTitleY("E' [GeV]");

H2F Q2_vs_xB = new H2F("Q2_vs_xB", "Q2 vs xB", 500, 0, 1, 500, 0, 13);
Q2_vs_xB.setTitleX("xB");
Q2_vs_xB.setTitleY("Q^2 [GeV^2]");

H2F W_vs_xB = new H2F("W_vs_xB", "W vs xB", 500, 0, 0.81, 500, wmin, wmax+0.5);
W_vs_xB.setTitleX("xB");
W_vs_xB.setTitleY("W [GeV]");

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 500, wmin, wmax, 500, -phimax, phimax);
Phi_vs_W.setTitleX("W [GeV]");
Phi_vs_W.setTitleY("Phi [deg]");

H2F xsect_vs_xB = new H2F("xsect_vs_xB", "generating #sigma vs xB", 500, 0.0, 0.81, 500, -1, 150);
xsect_vs_xB.setTitleX("xB");
xsect_vs_xB.setTitleY("#sigma");


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
    int pid_gen = 0;
    byte q_gen = 0;
    float px_gen = 0;
    float py_gen = 0;
    float pz_gen = 0;
    float beta_gen =0;
    float mom_gen = 0;
    double phi_gen = 0;
    double theta_gen =  0;
    float vz_gen = 0;
    double Q2_gen = 0; 
    double W_gen =0;
    double E_prime_gen = 0;
    double xB_gen = 0;
    
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
        
        // get generated data
       if ( event.hasBank("MC::Particle") ) {
            DataBank bank_gen = event.getBank("MC::Particle");
            
            for (int k = 0; k < bank_gen.rows(); k++) {
                // get values
                px_gen = bank_gen.getFloat("px", k);
                py_gen = bank_gen.getFloat("py", k);
                pz_gen = bank_gen.getFloat("pz", k);
                
                // calculate values
                mom_gen = (float) Math.sqrt(px_gen * px_gen + py_gen * py_gen + pz_gen * pz_gen);
                phi_gen = Math.atan2((double) py_gen,(double) px_gen);
                theta_gen = Math.acos((double) pz_gen/(double) mom_gen);
                
                Vector3 e_vec_3_gen = new Vector3(px_gen, py_gen, pz_gen); //3 vector e'
                LorentzVector e_vec_prime_gen = new LorentzVector(); //4 vector e'
                e_vec_prime_gen.setVectM(e_vec_3_gen, e_mass);
                
                LorentzVector q_vec_gen = new LorentzVector(); //4 vector q
                q_vec_gen.copy(e_vec); //e - e'
                q_vec_gen.sub(e_vec_prime_gen);
                Q2_gen = -q_vec_gen.mass2(); //-q^2
                
                LorentzVector w_vec_gen = new LorentzVector(); //4 vector used to calculate W
                w_vec_gen.copy(p_vec); //p + q
                w_vec_gen.add(q_vec_gen);
                W_gen = w_vec_gen.mass();
                
                E_prime_gen = e_vec_prime_gen.e();
                xB_gen = Q2_gen/(2.0*p_mass*(en-E_prime_gen));
                
                theta_gen *= 180/Math.PI;
                phi_gen *= 180/Math.PI;
                
                // Fill histos
                theta_hist_gen.fill(theta_gen);
                phi_hist_gen.fill(phi_gen);
                mom_hist_gen.fill(mom_gen);
                
                W_hist_gen.fill(W_gen);
                Q2_hist_gen.fill(Q2_gen);
                xB_hist_gen.fill(xB_gen);
            }
            
        }
        
        // get reconstructed data
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj") && event.hasBank("MC::Event")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            DataBank bank_mcEvent = event.getBank("MC::Event");
            
            DataBank bank_gen = event.getBank("MC::Particle");
            
            //System.out.println("number of bank_rec rows: " + bank_rec.rows() + ", # of gen bank rows: " + bank_gen.rows() );
            
            for (int k = 0; k < bank_rec.rows(); k++) tot_xsect += bank_mcEvent.getFloat("weight", 0);
                
            for (int k = 0; k < bank_rec.rows(); k++) {
   
                weight = bank_mcEvent.getFloat("weight", 0);
    
                pid = bank_rec.getInt("pid", k);
                q = bank_rec.getByte("charge", k);
                px = bank_rec.getFloat("px", k);
                py = bank_rec.getFloat("py", k);
                pz = bank_rec.getFloat("pz", k);
                beta = bank_rec.getFloat("beta", k);
    
                // get gen values
                px_gen = bank_gen.getFloat("px", 0);
                py_gen = bank_gen.getFloat("py", 0);
                pz_gen = bank_gen.getFloat("pz", 0);
                
                // calculate values
                mom_gen = (float) Math.sqrt(px_gen * px_gen + py_gen * py_gen + pz_gen * pz_gen);
                phi_gen = Math.atan2((double) py_gen,(double) px_gen);
                theta_gen = Math.acos((double) pz_gen/(double) mom_gen);
                
                
                mom = (float) Math.sqrt(px * px + py * py + pz * pz);
                phi = Math.atan2((double) py,(double) px);
                theta = Math.acos((double) pz/(double) mom);
    
                theta_gen *= 180/Math.PI;
                phi_gen *= 180/Math.PI;
                
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                vz = bank_rec.getFloat("vz", k);
    
                // pick electrons
                //if (pid != 11) continue;
                if(q != -1) continue;
                
                Vector3 e_vec_3_gen = new Vector3(px_gen, py_gen, pz_gen); //3 vector e'
                LorentzVector e_vec_prime_gen = new LorentzVector(); //4 vector e'
                e_vec_prime_gen.setVectM(e_vec_3_gen, e_mass);
                
                LorentzVector q_vec_gen = new LorentzVector(); //4 vector q
                q_vec_gen.copy(e_vec); //e - e'
                q_vec_gen.sub(e_vec_prime_gen);
                Q2_gen = -q_vec_gen.mass2(); //-q^2
                
                LorentzVector w_vec_gen = new LorentzVector(); //4 vector used to calculate W
                w_vec_gen.copy(p_vec); //p + q
                w_vec_gen.add(q_vec_gen);
                W_gen = w_vec_gen.mass();
                
                E_prime_gen = e_vec_prime_gen.e();
                xB_gen = Q2_gen/(2.0*p_mass*(en-E_prime_gen));
                
    
                Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
                LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
                e_vec_prime.setVectM(e_vec_3, e_mass);
    
                // ---------------- Cut on E' and theta -------------------
                //if(e_vec_prime.e() < 0.1 * en){continue;}  //cut below 10% beam
                //if(theta < 5 || theta > 40){continue;}     //cut outside of 5 and 40 degrees for FD
                // --------------------------------------------------------

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
                
                // ------------------------ Cuts --------------------------
                //if(W < 2) continue;                    // cut below 2 GeV/c^2
                
                // --------------------------------------------------------
                
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
                
                Q2_vs_W.fill(W,Q2);
                Phi_vs_W.fill(W,phi);
                E_vs_Theta.fill(theta,e_vec_prime.e()); // check this with calculated E'
                Q2_vs_xB.fill(xB,Q2);
                W_vs_xB.fill(xB,W);
                //if (Q2>1 && Q2 <3) 
                xsect_vs_xB.fill(xB,weight);
    
                // Calculate max values of each param
                if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} 
                if(theta > thetamax){thetamax = theta;}
                if(phi > phimax){phimax = phi;}
                if(vz > vzmax){vzmax = vz;}
                
                theta_hist_res.fill(theta_gen-theta);
                phi_hist_res.fill(phi_gen-phi);
                mom_hist_res.fill(mom_gen-mom);
                
                W_hist_res.fill(W_gen-W);
                Q2_hist_res.fill(Q2_gen-Q2);
                xB_hist_res.fill(xB_gen-xB);

            } // end for
        } // end if
    } // end while
    
    
    //mom_hist_res.sub(momentum);
    
} // end open file

TCanvas can_1d = new TCanvas("can", 1100, 600);
can_1d.setTitle("Theta, phi & mom resolutions");
can_1d.divide(3,3);
can_1d.cd(0);
can_1d.draw(theta_hist_gen);
can_1d.cd(1);
can_1d.draw(phi_hist_gen);
can_1d.cd(2);
can_1d.draw(mom_hist_gen);
can_1d.cd(3);
can_1d.draw(theta_hist); 
can_1d.cd(4);
can_1d.draw(phi_hist);
can_1d.cd(5);
can_1d.draw(momentum);
can_1d.cd(6);
can_1d.draw(theta_hist_res);
can_1d.cd(7);
can_1d.draw(phi_hist_res);
can_1d.cd(8);
can_1d.draw(mom_hist_res);

//for(int k =0; k < 8; k++){ 
//    xB_histmap.get(k).setLineColor(k);
//    can_1d.draw(xB_histmap.get(k),"same");
    //can.getPad().setLegend(true);
    //can.getPad().setLegendPosition(20, 20);
    //legend.AddEntry(xB_histmap.get(k),"Q2 < " + 1*k+1,"l");
//}
can_1d.save("figs/uncut/1D_spectra.png");

TCanvas can_1d_b = new TCanvas("can", 1100, 600);
can_1d_b.setTitle("W, Q2, xB resolutions");
can_1d_b.divide(3,3);
can_1d_b.cd(0);
can_1d_b.draw(W_hist_gen);
can_1d_b.cd(1);
can_1d_b.draw(Q2_hist_gen);
can_1d_b.cd(2);
can_1d_b.draw(xB_hist_gen);
can_1d_b.cd(3);
can_1d_b.draw(W_hist); 
can_1d_b.cd(4);
can_1d_b.draw(Q2_hist);
can_1d_b.cd(5);
can_1d_b.draw(xB_hist);
can_1d_b.cd(6);
can_1d_b.draw(W_hist_res);
can_1d_b.cd(7);
can_1d_b.draw(Q2_hist_res);
can_1d_b.cd(8);
can_1d_b.draw(xB_hist_res);
can_1d_b.save("figs/uncut/1D_kin_spectra.png");

TCanvas can_2d = new TCanvas("can", 1100, 600);
can_2d.divide(2,2);
can_2d.cd(0);
can_2d.draw(Q2_vs_W);
can_2d.cd(1);
can_2d.draw(E_vs_Theta);
can_2d.cd(2);
can_2d.draw(Q2_vs_xB);
can_2d.cd(3);
can_2d.draw(W_vs_xB);
can_2d.save("figs/uncut/2d_spectra.png");

/*
HashMap<Integer,TCanvas> canvasmap = new HashMap<Integer,TCanvas>();
for(int i : histmap.keySet()){
canvasmap.put(i, new TCanvas("can", 800,600));
canvasmap.get(i).draw(histmap.get(i));
canvasmap.get(i).save("phivstheta" + i + ".png");
}*/