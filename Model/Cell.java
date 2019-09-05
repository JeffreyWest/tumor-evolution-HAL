package Model;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;

public class Cell extends AgentSQ2Dunstackable<TumorEvolution> {

    Clone clone;

    Cell Init(Clone parent){
        this.clone=parent;
        this.clone.IncPop(); // always increment after init.
        return this;
    }


    Cell Mutate(){

        // if driver mutation, add one to # of drivers
        int new_kd = (G.rn.Double() < ( G.mu_d)) ? clone.kd+1 : clone.kd;

        // if passenger mutation, add one to # of passengers
        int new_kp = (G.rn.Double() <(G.mu_p)) ? clone.kp+1 : clone.kp;

        if ((new_kp != clone.kp) || (new_kd != clone.kd)) {
            // mutation has occured
            this.clone.DecPop();
            this.clone = new Clone(this.clone,new_kd,new_kp);      // update clonal lineage
            this.clone.IncPop();
        }

        return this;
    }

    void Birth(){
        // Passengers lower birth rate; Drivers raise birth rate
        double clone_specific_birth_rate = Math.pow(1.0+G.sd,(double)clone.kd)*G.birth_rate;

        if(G.rn.Double()<(clone_specific_birth_rate)) {
            int nDivOptions = G.MapEmptyHood(G.neighborhood, Xsq(), Ysq());
            if (nDivOptions == 0) {
                return; // no space to divide
            }

            int nextAgentID = G.neighborhood[G.rn.Int(nDivOptions)];

            // create daughter
            Cell daughter = G.NewAgentSQ(nextAgentID).Init(clone);//, this.progenyID, this.parentID);

            daughter.Mutate();
        }

        return;

    }

    // constant death rate
    void Death() {
        if(G.rn.Double()<(G.death_rate )){
            this.clone.DecPop();
            Dispose();
        }
    }


}
