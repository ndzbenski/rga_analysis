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
double wmin = 1.8;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

HipoDataSource reader = new HipoDataSource();

// The 1D histos
H1F theta_hist = new H1F("theta", "theta", 500, 0, thetamax+1);
theta_hist.setTitleX("theta [deg]");

H1F phi_hist = new H1F("phi", "phi", 500, 0, 370);
phi_hist.setTitleX("phi [deg]");

H1F momentum = new H1F("momentum", "momentum", 500, 0, 11);
momentum.setTitleX("momentum [GeV]");

H1F W_hist = new H1F("W", "W", 500, 0, wmax);
W_hist.setTitleX("W [GeV]");

H1F Q2_hist = new H1F("Q2", "Q2", 50, 0, 13);
Q2_hist.setTitleX("Q^2 [GeV^2]");

H1F Eprime_hist = new H1F("Eprime", "E'", 50, 0, 13);
Eprime_hist.setTitleX("E' [GeV]");

// H1F xB_hist = new H1F("xB", "xB", 100, 0, 0.9);
//xB_hist.setTitleX("xB");

// Create multiple histograms for each Q2 bin 
HashMap<Integer,H1F> xB_histmap = new HashMap<Integer,H1F>();
for(int i = 0; i < 8; i++){xB_histmap.put(i,new H1F("xB", 100, 0, 0.9));}

H1F xsection_hist = new H1F("xsection", "Generating Cross Section (#sigma)", 100, 0, 1);
xsection_hist.setTitleX("#sigma_{gen}");

// 2D Histos
H2F Q2_vs_W = new H2F("Q2_vs_W", "Q2 vs W", 500, wmin, wmax, 500, 0.0, 13);
Q2_vs_W.setTitleX("W [GeV]");
Q2_vs_W.setTitleY("Q^2 [GeV^2]");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E' vs Theta", 500, 5, thetamax+1, 500, 0, enmax);
E_vs_Theta.setTitleX("Theta [deg]");
E_vs_Theta.setTitleY("E' [GeV]");

H2F Q2_vs_xB = new H2F("Q2_vs_xB", "Q2 vs xB", 500, 0, 1, 500, 0, 13);
Q2_vs_xB.setTitleX("xB");
Q2_vs_xB.setTitleY("Q^2 [GeV^2]");

H2F W_vs_xB = new H2F("W_vs_xB", "W vs xB", 500, 0, 0.81, 500, wmin, wmax);
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
    
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        double tot_xsect = 0;
        
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj") && event.hasBank("MC::Event")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            DataBank bank_mcEvent = event.getBank("MC::Event");
            
            for (int k = 0; k < bank_rec.rows(); k++) tot_xsect += bank_mcEvent.getFloat("weight", 0);
                
            for (int k = 0; k < bank_rec.rows(); k++) {
   
                float weight = bank_mcEvent.getFloat("weight", 0);
    
                int pid = bank_rec.getInt("pid", k);
                byte q = bank_rec.getByte("charge", k);
                float px = bank_rec.getFloat("px", k);
                float py = bank_rec.getFloat("py", k);
                float pz = bank_rec.getFloat("pz", k);
                float beta = bank_rec.getFloat("beta", k);
    
                float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
                double phi = Math.atan2((double) py,(double) px);
                double theta = Math.acos((double) pz/(double) mom);
    
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                float vz = bank_rec.getFloat("vz", k);
    
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

                momentum.fill(mom);
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                double Q2 = -q_vec.mass2(); //-q^2
                
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                double W = w_vec.mass();
                //double W = p_mass*p_mass + 2.0*p_mass*(en-E_prime) -Q2
    
                //double E_prime = Q2/(4.0*en*(Math.sin(theta*Math.PI/360.0)));
                double E_prime = e_vec_prime.e();
                double xB = Q2/(2.0*p_mass*(en-E_prime));
                
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

            } // end for
        } // end if
    } // end while
} // end open file

TCanvas can_1d = new TCanvas("can", 1100, 600);
can_1d.title("Reconstructed");
can_1d.divide(3,3);
can_1d.cd(1);
can_1d.draw(theta_hist);
can_1d.cd(2);
can_1d.draw(phi_hist);
can_1d.cd(3);
can_1d.draw(momentum);
can_1d.cd(4);
can_1d.draw(W_hist); 
can_1d.cd(5);
can_1d.draw(Q2_hist);
can_1d.cd(6);
for(int k =0; k < 8; k++){ 
    xB_histmap.get(k).setLineColor(k);
    can_1d.draw(xB_histmap.get(k),"same");
    //can.getPad().setLegend(true);
    //can.getPad().setLegendPosition(20, 20);
    //legend.AddEntry(xB_histmap.get(k),"Q2 < " + 1*k+1,"l");
}
can_1d.save("figs/1D_spectra.png");

//TCanvas can_ = new TCanvas("can", 1100, 600);

TCanvas can_2d = new TCanvas("can", 1100, 600);
can_2d.divide(3,2);
can_2d.cd(1);
can_2d.draw(Q2_vs_W);
can_2d.cd(2);
can_2d.draw(E_vs_Theta);
can_2d.cd(3);
can_2d.draw(Q2_vs_xB);
can_2d.cd(4);
can_2d.draw(W_vs_xB);
can_2d.cd(5);
can_2d.draw(xsect_vs_xB);
can_2d.save("figs/2d_spectra.png");

/*
HashMap<Integer,TCanvas> canvasmap = new HashMap<Integer,TCanvas>();
for(int i : histmap.keySet()){
canvasmap.put(i, new TCanvas("can", 800,600));
canvasmap.get(i).draw(histmap.get(i));
canvasmap.get(i).save("phivstheta" + i + ".png");
}*/