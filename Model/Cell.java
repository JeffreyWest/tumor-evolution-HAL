package Model;

import Framework.GridsAndAgents.AgentSQ2Dunstackable;

public class Cell extends AgentSQ2Dunstackable<TumorEvolution> {
    int kd;
    int kp;
    int progenyID;
    int parentID;

    Cell Init(int kp0, int kd0, int progenyID0, int parentID0){
        kp = kp0;
        kd = kd0;
        progenyID = progenyID0;
        parentID = parentID0;
        G.driver_status[0] = 1;

        return this;
    }

    Cell Mutate(){
        boolean mutated = false;
        boolean driver_mutated = false;

        // driver mutation
        if((G.rn.Double() < ( G.mu_d))) {
            kd++;
            if (kd > G.KdMAX) { G.KdMAX++; }
            mutated = true;
            driver_mutated = true;
        }

        // passenger mutation
        if((G.rn.Double() <(G.mu_p))) {
            kp++;
            if (kp > G.KpMAX) { G.KpMAX++; }
            mutated = true;
        }

        if (mutated) {
            parentID = progenyID;
            progenyID = G.progenyNextID;
            G.progenyToParentIDs[progenyID] = parentID;
            G.driver_status[progenyID] = (driver_mutated) ? kd : kd;
            G.passenger_number[progenyID] = kp;
            G.progenyNextID++;

            if (G.progenyNextID >= G.ExpectedNumberOfClones) {
                System.out.println("\nERROR :::: increase your expected number of clones: " + Integer.toString(G.ExpectedNumberOfClones) + " is too small.\nSet this is in the TumorEvolution class.");
                System.exit(1);
            }

        }


        return this;
    }

    Cell Birth(){
        // Passengers lower birth rate; Drivers raise birth rate
        double effective_birth_rate = Math.pow(1.0+G.sd,(double)kd)*G.birth_rate;

        if(G.rn.Double()<(effective_birth_rate)) {
            int nDivOptions = G.MapEmptyHood(G.neighborhood, Xsq(), Ysq());
            if (nDivOptions == 0) {
                return null; // no space to divide
            }
            int nextAgentID = G.neighborhood[G.rn.Int(nDivOptions)];
            return G.NewAgentSQ(nextAgentID).Init(this.kp, this.kd, this.progenyID, this.parentID).Mutate();

        }else {
            return null;
        }
    }

    // constant death rate
    void Death() {
        if(G.rn.Double()<(G.death_rate )){
            Dispose();
        }
    }
}
