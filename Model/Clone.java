package Model;
import HAL.Tools.PhylogenyTracker.Genome;

import static HAL.Util.*;

class Clone extends Genome<Clone> {

    // clone-specific attributes
    int color;
    int kd;
    int kp;

    public Clone(Clone parent, int kd, int kp) {
        super(parent, false);

        // set clone-specific attributes
        this.kd = kd;
        this.kp = kp;
        this.color = CategorialColor((kd - 1) % 19 );
    }

    public String Hex() {
        return String.format("#%02x%02x%02x", GetRed256(color), GetGreen256(color), GetBlue256(color));
    }


}