import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// usage:
// rungroovy uncut_analysis.groovy <mc_data_file.hipo> <rga_data_file.hipo> <energy>
// for use on RGA data

double en = Double.parseDouble(args[2]);
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

// The MC reconstructed 1D histos
H1F theta_hist_mc = new H1F("theta_mc", "theta_mc", 100, 0, thetamax+5);
theta_hist_mc.setTitleX("theta_mc [deg]");

H1F phi_hist_mc = new H1F("phi_mc", "phi_mc", 100, -phimax, phimax);
phi_hist_mc.setTitleX("phi_mc [deg]");

H1F mom_hist_mc = new H1F("momentum mc", "momentum mc", 100, 0, 11);
mom_hist_mc.setTitleX("momentum_mc [GeV]");

H1F W_hist_mc = new H1F("W_mc", "W_mc", 100, wmin, wmax+0.5);
W_hist_mc.setTitleX("W_mc [GeV]");

H1F Q2_hist_mc = new H1F("Q2_mc", "Q2_mc", 50, 0, 13);
Q2_hist_mc.setTitleX("Q^2 mc [GeV^2]");

H1F xB_hist_mc = new H1F("xB_mc", "xB_mc", 100, 0, 1);
xB_hist_mc.setTitleX("xB_mc");


// The reconstructed 1D histos 
H1F theta_hist = new H1F("theta", "theta", 100, 0, thetamax+5);
theta_hist.setTitleX("theta [deg]");

H1F phi_hist = new H1F("phi", "phi", 100, -phimax, phimax);
phi_hist.setTitleX("phi [deg]");

H1F momentum = new H1F("momentum", "momentum", 100, 0, 11);
momentum.setTitleX("momentum [GeV]");

H1F W_hist = new H1F("W", "W", 100, wmin, wmax+0.5);
W_hist.setTitleX("W [GeV]");

H1F Q2_hist = new H1F("Q2", "Q2", 50, 0, 13);
Q2_hist.setTitleX("Q^2 [GeV^2]");

H1F Eprime_hist = new H1F("Eprime", "E'", 50, 0, 13);
Eprime_hist.setTitleX("E' [GeV]");

H1F xB_hist = new H1F("xB", "xB", 100, 0, 1);
xB_hist.setTitleX("xB");


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

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 100, wmin, wmax, 100, -phimax, phimax);
Phi_vs_W.setTitleX("W [GeV]");
Phi_vs_W.setTitleY("Phi [deg]");

H2F xsect_vs_xB = new H2F("xsect_vs_xB", "generating #sigma vs xB", 100, 0.0, 0.81, 100, -1, 150);
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

// for MC data
// reconstructed data variables
new File('.', args[0]).eachLine { line ->
    reader.open(line);
    
    float weight = 0;
    int pid_mc = 0;
    float px_mc = 0;
    float py_mc = 0;
    float pz_mc = 0;
    float mom_mc = 0;
    double phi_mc = 0;
    double theta_mc =  0;
    double Q2_mc = 0; 
    double W_mc =0;
    double E_prime_mc = 0;
    double xB_mc = 0;
    
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
         // get MC reconstructed data
       if ( event.hasBank("RECHB::Particle") ) {
            DataBank bank_mc = event.getBank("RECHB::Particle");
            
            for (int k = 0; k < bank_mc.rows(); k++) {
                // get values
                px_mc = bank_mc.getFloat("px", k);
                py_mc = bank_mc.getFloat("py", k);
                pz_mc = bank_mc.getFloat("pz", k);
                pid_mc = bank_mc.getInt("pid", k);
                
                if(pid_mc != 11) continue;
                
                // calculate values
                mom_mc = (float) Math.sqrt(px_mc * px_mc + py_mc * py_mc + pz_mc * pz_mc);
                phi_mc = Math.atan2((double) py_mc,(double) px_mc);
                theta_mc = Math.acos((double) pz_mc/(double) mom_mc);
                
                Vector3 e_vec_3_mc = new Vector3(px_mc, py_mc, pz_mc); //3 vector e'
                LorentzVector e_vec_prime_mc = new LorentzVector(); //4 vector e'
                e_vec_prime_mc.setVectM(e_vec_3_mc, e_mass);
                
                LorentzVector q_vec_mc = new LorentzVector(); //4 vector q
                q_vec_mc.copy(e_vec); //e - e'
                q_vec_mc.sub(e_vec_prime_mc);
                Q2_mc = -q_vec_mc.mass2(); //-q^2
                
                LorentzVector w_vec_mc = new LorentzVector(); //4 vector used to calculate W
                w_vec_mc.copy(p_vec); //p + q
                w_vec_mc.add(q_vec_mc);
                W_mc = w_vec_mc.mass();
                
                E_prime_mc = e_vec_prime_mc.e();
                xB_mc = Q2_mc/(2.0*p_mass*(en-E_prime_mc));
                
                theta_mc *= 180/Math.PI;
                phi_mc *= 180/Math.PI;
                
                // Fill histos
                theta_hist_mc.fill(theta_mc);
                theta_hist_mc.setLineColor(3); 
                phi_hist_mc.fill(phi_mc);
                phi_hist_mc.setLineColor(3); 
                mom_hist_mc.fill(mom_mc);
                mom_hist_mc.setLineColor(3); 
                
                W_hist_mc.fill(W_mc);
                W_hist_mc.setLineColor(3); 
                Q2_hist_mc.fill(Q2_mc);
                Q2_hist_mc.setLineColor(3); 
                xB_hist_mc.fill(xB_mc);
                xB_hist_mc.setLineColor(3); 
                
            } // end for loop
        }// end if 
    } // end event
} // end open MC data file
    
new File('.', args[1]).eachLine { line ->
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
        
        // get reconstructed data
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            
            //System.out.println("number of bank_rec rows: " + bank_rec.rows() + ", # of gen bank rows: " + bank_mc.rows() );
                
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
    
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                Q2 = -q_vec.mass2(); //-q^2
                
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                W = w_vec.mass();
                
                E_prime = e_vec_prime.e();
                xB = Q2/(2.0*p_mass*(en-E_prime));
                
                // Fill histos
                theta_hist.fill(theta);
                phi_hist.fill(phi);
                momentum.fill(mom);
                
                xB_hist.fill(xB);
                Q2_hist.fill(Q2);
                W_hist.fill(W);
                
                Eprime_hist.fill(E_prime);
                
                Q2_vs_W.fill(W,Q2);
                Phi_vs_W.fill(W,phi);
                E_vs_Theta.fill(theta,e_vec_prime.e());
                Q2_vs_xB.fill(xB,Q2);
                W_vs_xB.fill(xB,W);
    
                // Calculate max values of each param
                if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} 
                if(theta > thetamax){thetamax = theta;}
                if(phi > phimax){phimax = phi;}
                if(vz > vzmax){vzmax = vz;}
                
            } // end for
        } // end if
    } // end while
    
    
    //mom_hist_res.sub(momentum);
    
} // end open file

TCanvas can_1d = new TCanvas("can", 1100, 600);
can_1d.setTitle("Theta, phi & mom resolutions");
can_1d.divide(3,2);
can_1d.draw(theta_hist);
can_1d.draw(theta_hist_mc, "same");
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20);
can_1d.cd(1);
can_1d.draw(phi_hist);
can_1d.draw(phi_hist_mc, "same");
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20); 
can_1d.cd(2);
can_1d.draw(momentum);
can_1d.draw(mom_hist_mc, "same");
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20); 
can_1d.cd(3)
can_1d.draw(W_hist); 
can_1d.draw(W_hist_mc,"same"); 
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20);
can_1d.cd(4);
can_1d.draw(Q2_hist);
can_1d.draw(Q2_hist_mc,"same"); 
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20);
can_1d.cd(5);
can_1d.draw(xB_hist);
can_1d.draw(xB_hist_mc,"same"); 
can_1d.getPad().setLegend(true);
can_1d.getPad().setLegendPosition(20, 20);
can_1d.save("figs/uncut/1D_spectra.png");

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