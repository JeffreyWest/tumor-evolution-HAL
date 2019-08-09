package Framework.GridsAndAgents;

import java.io.Serializable;

/**
 * AgentBase functions are shared amongst all agent types
 */
public abstract class AgentBase<T> implements Serializable{
    int stateID;
    int iList;
    int birthTick;
    boolean alive;
    /**
     * returns the grid that the agent belongs to (this is a permanent agent property, not a function)
     */
    public final T G;
    AgentListNode myNodes;

    public AgentBase(){
        G=null;
    }

    /**
     * Returns whether the agent is alive or has been disposed
     */
    public boolean IsAlive(){
        return alive;
    }

    /**
     * returns the tick on which the agent was born
     */

    public int BirthTick(){
        return birthTick;
    }

    /**
     * returns the age of the agent, in ticks. Be sure to use IncTick on the AgentGrid appropriately for this function to work.
     */
    abstract public int Age();

    /**
     * removes the agent from the grid
     */
    abstract public void Dispose();

}
