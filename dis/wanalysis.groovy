import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// usage:
// rungroovy wanalysis.groovy <filename.hipo> <energy>

double en = Double.parseDouble(args[1]);
double enmax = en+0.1; //GeV
double thetamax = 40;  //degrees
double phimax = 180;   //degrees
double vzmax = 50;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

HipoDataSource reader = new HipoDataSource();

H1F momentum = new H1F("momentum", "momentum", 500, 0, 10);
momentum.setTitleX("momentum");

H1F W_hist = new H1F("W", "W", 500, 0, wmax);
W_hist.setTitleX("W");

H1F Q2_hist = new H1F("Q2", "Q2", 50, 0, 10);
Q2_hist.setTitleX("Q2");

H1F xB_hist = new H1F("xB", "xB", 50, 0, 1);
xB_hist.setTitleX("xB");

H2F W_vs_Q2 = new H2F("W_vs_Q2", "W_vs_Q2", 500, 0.0, enmax, 500, 0.0, wmax);
W_vs_Q2.setTitleX("Q2");
W_vs_Q2.setTitleY("W");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E_vs_Theta", 500, 5, thetamax, 500, 0, enmax);
E_vs_Theta.setTitleX("Theta");
E_vs_Theta.setTitleY("E'");

H2F Q2_vs_xB = new H2F("Q2_vs_xB", "Q2_vs_xB", 500, 0, 1, 500, 0, 15);
Q2_vs_xB.setTitleX("xB");
Q2_vs_xB.setTitleY("Q2");


H2F W_vs_xB = new H2F("W_vs_xB", "W_vs_xB", 500, 0, 1, 500, 0, wmax);
W_vs_xB.setTitleX("xB");
W_vs_xB.setTitleY("W");

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 500, 0, wmax, 500, -phimax, phimax);
Phi_vs_W.setTitleX("W");
Phi_vs_W.setTitleY("Phi");


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
    
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj") && event.hasBank("MC::Event")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            DataBank bank_mcEvent = event.getBank("MC::Event");
            
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
    
                if(e_vec_prime.e() < 0.1 * en){continue;} //cut below 10% beam
                if(theta < 5 || theta > 40){continue;} //cut outside of 5 and 40 degrees for FD

                // Fill Histograms
                momentum.fill(mom);
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                double Q2 = -q_vec.mass2(); //-q^2
                
    
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                double W = w_vec.mass();
    
                // xB = Q^2/2M_p(E-E')
                double E_prime = Q2/(4.0*en*(Math.sin(theta*Math.PI/360.0)));
                double xB = Q2/(2.0*p_mass*(en-E_prime));
                
                // Make cut on W < 2 MeV/c^2
                if(W < 2) continue;
    
                //double W = p_mass*p_mass + 2.0*p_mass*(en-E_prime) -Q2
                Q2_hist.fill(Q2);
                xB_hist.fill(xB);
                W_hist.fill(W);
                W_vs_Q2.fill(Q2,W);
                Phi_vs_W.fill(W,phi);
    
                if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} //calculate max values of each param
                if(theta > thetamax){thetamax = theta;}
                if(phi > phimax){phimax = phi;}
                if(vz > vzmax){vzmax = vz;}
    
                E_vs_Theta.fill(theta,e_vec_prime.e());
                Q2_vs_xB.fill(xB,Q2);
                W_vs_xB.fill(xB,W);
    
    
            } // end for
        } // end if
    } // end while
} // end open file


TCanvas can2 = new TCanvas("can", 800, 600);
can2.draw(Q2_hist);
can2.save("figs/Q2.png");

TCanvas can3 = new TCanvas("can", 800, 600);
can3.draw(W_hist);
can3.save("figs/W.png");

TCanvas can8 = new TCanvas("can", 800, 600);
can8.draw(xB_hist);
can8.save("figs/xB.png");


TCanvas can = new TCanvas("can", 800, 600);
can.draw(W_vs_Q2);
can.save("figs/W_vs_Q2.png");


TCanvas can1 = new TCanvas("can", 800, 600);
can1.draw(momentum);
can1.save("figs/mom.png");

TCanvas can4 = new TCanvas("can", 800, 600);
can4.draw(E_vs_Theta);
can4.save("figs/EvsTheta.png");

TCanvas can5 = new TCanvas("can", 800, 600);
can5.draw(Q2_vs_xB);
can5.save("figs/Q2_vs_xB.png");


TCanvas can6 = new TCanvas("can", 800, 600);
can6.draw(W_vs_xB);
can6.save("figs/W_vs_xB.png");

/*
TCanvas can7 = new TCanvas("can", 800, 600);
can7.draw(Phi_vs_W);
can7.save("figs/PhivsW.png"); */



/*
HashMap<Integer,TCanvas> canvasmap = new HashMap<Integer,TCanvas>();
for(int i : histmap.keySet()){
canvasmap.put(i, new TCanvas("can", 800,600));
canvasmap.get(i).draw(histmap.get(i));
canvasmap.get(i).save("phivstheta" + i + ".png");
}*/