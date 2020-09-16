
function matrix_invert(M){
    // I use Guassian Elimination to calculate the inverse:
    // (1) 'augment' the matrix (left) by the identity (on the right)
    // (2) Turn the matrix on the left into the identity by elemetry row ops
    // (3) The matrix on the right is the inverse (was the identity matrix)
    // There are 3 elemtary row ops: (I combine b and c in my code)
    // (a) Swap 2 rows
    // (b) Multiply a row by a scalar
    // (c) Add 2 rows
    
    //if the matrix isn't square: exit (error)
    if(M.length !== M[0].length){return;}
    
    //create the identity matrix (I), and a copy (C) of the original
    var i=0, ii=0, j=0, dim=M.length, e=0, t=0;
    var I = [], C = [];
    for(i=0; i<dim; i+=1){
        // Create the row
        I[I.length]=[];
        C[C.length]=[];
        for(j=0; j<dim; j+=1){
            
            //if we're on the diagonal, put a 1 (for identity)
            if(i==j){ I[i][j] = 1; }
            else{ I[i][j] = 0; }
            
            // Also, make the copy of the original
            C[i][j] = M[i][j];
        }
    }
    
    // Perform elementary row operations
    for(i=0; i<dim; i+=1){
        // get the element e on the diagonal
        e = C[i][i];
        
        // if we have a 0 on the diagonal (we'll need to swap with a lower row)
        if(e==0){
            //look through every row below the i'th row
            for(ii=i+1; ii<dim; ii+=1){
                //if the ii'th row has a non-0 in the i'th col
                if(C[ii][i] != 0){
                    //it would make the diagonal have a non-0 so swap it
                    for(j=0; j<dim; j++){
                        e = C[i][j];       //temp store i'th row
                        C[i][j] = C[ii][j];//replace i'th row by ii'th
                        C[ii][j] = e;      //repace ii'th by temp
                        e = I[i][j];       //temp store i'th row
                        I[i][j] = I[ii][j];//replace i'th row by ii'th
                        I[ii][j] = e;      //repace ii'th by temp
                    }
                    //don't bother checking other rows since we've swapped
                    break;
                }
            }
            //get the new diagonal
            e = C[i][i];
            //if it's still 0, not invertable (error)
            if(e==0){console.log('non-invertible'); return;}
        }
        
        // Scale this row down by e (so we have a 1 on the diagonal)
        for(j=0; j<dim; j++){
            C[i][j] = C[i][j]/e; //apply to original matrix
            I[i][j] = I[i][j]/e; //apply to identity
        }
        
        // Subtract this row (scaled appropriately for each row) from ALL of
        // the other rows so that there will be 0's in this column in the
        // rows above and below this one
        for(ii=0; ii<dim; ii++){
            // Only apply to other rows (we want a 1 on the diagonal)
            if(ii==i){continue;}
            
            // We want to change this element to 0
            e = C[ii][i];
            
            // Subtract (the row above(or below) scaled by e) from (the
            // current row) but start at the i'th column and assume all the
            // stuff left of diagonal is 0 (which it should be if we made this
            // algorithm correctly)
            for(j=0; j<dim; j++){
                C[ii][j] -= e*C[i][j]; //apply to original matrix
                I[ii][j] -= e*I[i][j]; //apply to identity
            }
        }
    }
    
    //we've done all operations, C should be the identity
    //matrix I should be the inverse:
    return I;
}

function matrix_multiply(a, b) {
    var aNumRows = a.length, aNumCols = a[0].length,
        bNumRows = b.length, bNumCols = b[0].length,
        m = new Array(aNumRows);  // initialize array of rows
    for (var r = 0; r < aNumRows; ++r) {
      m[r] = new Array(bNumCols); // initialize the current row
      for (var c = 0; c < bNumCols; ++c) {
        m[r][c] = 0;             // initialize the current cell
        for (var i = 0; i < aNumCols; ++i) {
          m[r][c] += a[r][i] * b[i][c];
        }
      }
    }
    return m;
}

function matrix_sum(a,b) {
    var result = [];
    result = new Array(a.length);
    for (var i = 0; i < result.length; i++) {
        result[i] = new Array(a[i].length);
        for (var j = 0; j < result[i].length; j++) {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

function matrix_scalar(s,m) {
    var result = [];
    result = new Array(m.length);
    for (var i = 0; i < result.length; i++) {
        result[i] = new Array(m[i].length);
        for (var j = 0; j < result[i].length; j++) {
            result[i][j] = s*m[i][j];
        }
    }
    return result;
}

function combine_arrays(x) {
    let res = [];
    for (var i = 0; i < x.length; i++) {
        for (var j = 0; j < x[i].length; j++) {
            res = res.concat(x[i][j]);
        }
    }
    return res;
}

function format_data(x,y,xf,yf) {
    if (!xf) {
        var xf = 1;
    }
    if (!yf) {
        var yf = 1;
    }
    let data = [];
    for (var i = 0; i < x.length; i++) {
        data = data.concat({x: xf*x[i], y: yf*y[i]});
    }
    return data;
}

function makeArr(a,b,n) {
    if(typeof n === "undefined") n = Math.max(Math.round(b-a)+1,1);
    if(n<2) { return n===1?[a]:[]; }
    var i,ret = Array(n);
    n--;
    for(i=n;i>=0;i--) { ret[i] = (i*b+(n-i)*a)/n; }
    return ret;
}

function quadratic_bezier(pts,Nr,type) {
    if (type == 'angles') {
        let N = pts[0];
        let A = pts[1];
        let E = pts[2];


        let m1 = Math.tan(A[0]);
        let m2 = Math.tan(A[1]);
        let c1 = N[1]-m1*N[0];
        let c2 = E[1]-m2*E[0];
        let Q = [(c2-c1)/(m1-m2),(m1*c2-m2*c1)/(m1-m2)];

        let t = makeArr(0,1,Nr);
        let x = matrix_scalar(1,t);
        let y = matrix_scalar(1,t);
        for (var i = 0; i < Nr; i++) {
            x[i] = Math.pow(1-t[i],2)*N[0] + 2*(1-t[i])*t[i]*Q[0] + t[i]*t[i]*E[0];
            y[i] = Math.pow(1-t[i],2)*N[1] + 2*(1-t[i])*t[i]*Q[1] + t[i]*t[i]*E[1];
        }
        var result = [x,y];
    }
    return result;
}

function AMRF(x,g) {
    let val = Math.pow((g+1)/2,-(g+1)/(2*g-2))*Math.pow(1+x*x*(g-1)/2,(g+1)/(2*g-2))/x;Math.pow((g+1)/2,-(g+1)/(2*g-2))*Math.pow(1+x*x*(g-1)/2,(g+1)/(2*g-2))/x;
    return val;
}

function AMR(x,g,sub) {
    let val = 0;
    if (!sub) {
        val = AMRF(x,g);
    } else if (sub == 'sub') {
        let M = 0.2;
        let dM = Math.pow(10,-6);
        let diff = 0;
        let F = AMRF(M,g) - x;
        while (Math.abs(F) > dM) {
            diff = (AMRF(M+dM,g)-AMRF(M,g))/dM;
            F = AMRF(M,g) - x;
            M = M - 0.01*F/diff;
        }
        val = M;
    } else if (sub == 'sup') {
        let M = 2;
        let dM = Math.pow(10,-6);
        let diff = 0;
        let F = AMRF(M,g) - x;
        while (Math.abs(F) > dM) {
            diff = (AMRF(M+dM,g)-AMRF(M,g))/dM;
            F = AMRF(M,g) - x;
            M = M - 0.01*F/diff;
        }
        val = M;
    }
    return val;
}

function try_shock(g,At,Ae,ANS) {
    let A1S = At;
    let M1 = AMR(ANS/A1S,g,'sup');
    let P1rP01 = Math.pow(1+M1*M1*(g-1)/2,-g/(g-1));
    let P2rP1 = 1 + 2*g/(g+1)*(M1*M1-1);
    let M2 = Math.sqrt((1+M1*M1*(g-1)/2)/(g*M1*M1-(g-1)/2));
    let P02rP2 = Math.pow(1+M2*M2*(g-1)/2,g/(g-1));
    let AErA2S = (Ae/At)*(A1S/ANS)*AMR(M2,g);
    let ME = AMR(AErA2S,g,'sub');
    let PErP0E = Math.pow(1+ME*ME*(g-1)/2,-g/(g-1));
    let PErP0 = PErP0E*P02rP2*P2rP1*P1rP01;
    return PErP0;
}


function find_state(P0,Pb,g,At,Ae) {

    /*Compute pressure ratios to find state*/
    let MEsub = AMR(Ae/At,g,'sub');
    let MEsup = AMR(Ae/At,g,'sup');
    let PErP0sub = Math.pow(1 + MEsub*MEsub*(g-1)/2,-g/(g-1));
    let PErP0sup = Math.pow(1 + MEsup*MEsup*(g-1)/2,-g/(g-1));
    let PErP0nse = (2*g)/(g+1)*(MEsup*MEsup-1)*PErP0sup;
    let PBrP0 = Pb/P0;
    let ANS = [];

    /* Check to find state */
    let state = [];
    if (PBrP0 >= PErP0sub) {
        state = 'subsonic';
        ANS = [MEsub, [PErP0sub,PErP0sub/PBrP0]];

    } else if (PBrP0 > PErP0nse) {
        state = 'shock';
        ANS = At*2;
        let dx = Math.pow(10,-6);
        let PErP0 = try_shock(g,At,Ae,ANS);
        let diff = 0;
        let ind = 0;

        while ((Math.abs(PErP0-PBrP0) > dx)&&(ind < 10000)) {
            PErP0 = try_shock(g,At,Ae,ANS);
            diff = (try_shock(g,At,Ae,ANS+dx)-try_shock(g,At,Ae,ANS))/dx;
            ANS = ANS - 0.1*(PErP0-PBrP0)/diff;
            ind = ind + 1;
        }

        let A1S = At;
        let M1 = AMR(ANS/A1S,g,'sup');
        let P1rP01 = Math.pow(1+M1*M1*(g-1)/2,-g/(g-1));
        let P2rP1 = 1 + 2*g/(g+1)*(M1*M1-1);
        let M2 = Math.sqrt((1+M1*M1*(g-1)/2)/(g*M1*M1-(g-1)/2));
        let AErA2S = (Ae/At)*(A1S/ANS)*AMR(M2,g);
        let ME = AMR(AErA2S,g,'sub');
        let P02rP01 = Math.pow(1+M2*M2*(g-1)/2,g/(g-1))*P2rP1*P1rP01;

        ANS = [ [M1,M2,ME], [P1rP01,P2rP1,PErP0,PErP0/PBrP0,P02rP01], [ANS/At,1/(At/Ae*AErA2S)]];

    } else if (PBrP0 > PErP0sup) {
        state = 'overexpanded';
        ANS = [MEsup,[PErP0sup,PErP0sup/PBrP0]];
    } else {
        state = 'underexpanded';
        ANS = [MEsup,[PErP0sup,PErP0sup/PBrP0]];
    }

    let results = [state,ANS];
    return results;
}

function bell_angle(E,Lf) {
    let ang1 = 18;
    let ang2 = 15;
    if ((E > 1.99)&&(E <= 100.05)&&(Lf <= 1.05)) {
        ang1 = (71-56/(0.01*E+1))+25*(1-Lf);
        ang2 = (23/E+3)+25*(1-Lf);
    } else if (E <= 1.99) {
        ang1 = 1.1*Math.atan(1/(Lf)*Math.tan(15*Math.PI/180))*180/Math.PI;
        ang2 = 0.9*Math.atan(1/(Lf)*Math.tan(15*Math.PI/180))*180/Math.PI;
    } else {
        E = 100;
        ang1 = (71-56/(0.01*E+1))+25*(1-Lf);
        ang2 = (23/E+3)+25*(1-Lf);
    }
    let ang = [ang1*Math.PI/180,ang2*Math.PI/180];
    return ang;
}

function nozzle_shape(type,At,E,Lf) {
    let x = [];
    let y = [];
    let Rt = Math.sqrt(At/Math.PI);
    let Re = Math.sqrt(E)*Rt;
    if (type == 'conical') {
        let Rc = 1.5*Rt;
        let Rs = 0.382*Rt;
        let x0c = -1*Rc*Math.sin(Math.PI/4);
        let a = 15*Math.PI/180;
        let x0s = Rs*Math.sin(a);
        let y0s = Rt+Rs*(1-Math.cos(a));
        let L = (Re-y0s)/Math.tan(a);
        let xc = makeArr(x0c,0,31);
        xc.pop();
        let xs = makeArr(0,x0s,11);
        xs.pop();
        let xd = makeArr(x0s,L+x0s,60);
        let yc = matrix_scalar(1,xc);
        let ind = 0;
        while (ind < xc.length) {
            yc[ind] = Rt + Math.sqrt(Rc*Rc + xc[ind]*xc[ind])-Rc;
            ind = ind + 1;
        }
        let ys = matrix_scalar(1,xs);
        ind = 0;
        while (ind < xs.length) {
            ys[ind] = Rt + Math.sqrt(Rs*Rs + xs[ind]*xs[ind])-Rs;
            ind = ind + 1;
        }
        let yd = matrix_scalar(1,xd);
        ind = 0;
        while (ind < xd.length) {
            yd[ind] = Math.tan(a)*(xd[ind]-x0s)+ys[ys.length-1];
            ind = ind + 1;
        }
        x = combine_arrays([xc,xs,xd]);
        y = combine_arrays([yc,ys,yd]);
    } else if (type == 'bell') {
        let Rc = 1.5*Rt;
        let Rs = 0.382*Rt;
        let x0c = -1*Rc*Math.sin(Math.PI/4);
        let [ai,af] = bell_angle(E,Lf);
        let x0s = Rs*Math.sin(ai);
        let y0s = Rt+Rs*(1-Math.cos(ai));
        let L = Lf*(Re-y0s)/Math.tan(15*Math.PI/180);
        let xc = makeArr(x0c,0,30);
        let xs = makeArr(0,x0s,11);
        xs.shift();
        let yc = matrix_scalar(1,xc);
        let ind = 0;
        while (ind < xc.length) {
            yc[ind] = Rt + Math.sqrt(Rc*Rc + xc[ind]*xc[ind])-Rc;
            ind = ind + 1;
        }
        let ys = matrix_scalar(1,xs);
        ind = 0;
        while (ind < xs.length) {
            ys[ind] = Rt + Math.sqrt(Rs*Rs + xs[ind]*xs[ind])-Rs;
            ind = ind + 1;
        }
        let [xd,yd] = quadratic_bezier([ [xs[xs.length-1], ys[ys.length-1]], [ai, af], [L+x0s,Re]],61,'angles');
        xd.shift();
        yd.shift();

        x = combine_arrays([xc,xs,xd]);
        y = combine_arrays([yc,ys,yd]);
    }

    let result = {x: x, y: y};
    return result;
}

function area_x(type,At,E,Lf) {
    let data = nozzle_shape(type,At,E,Lf);
    let x = data.x;
    let r = data.y;
    let area = matrix_scalar(1,r);
    for (var i = 0; i < r.length; i++) {
        area[i] = Math.PI*Math.pow(r[i],2);
    }
    let result = {x: x, a: area, r: r};
    return result;
}

function write_vectors(state_info,g,At,E,Lf,type,Pc,Tc,R) {
    let area = area_x(type,At,E,Lf);
    let M = matrix_scalar(1,area.a);
    let P = matrix_scalar(1,area.a);
    let T = matrix_scalar(1,area.a);
    let V = matrix_scalar(1,area.a);
    if ((state_info[0] == 'overexpanded')||(state_info[0] == 'underexpanded')) {
        let [state,[ME,[PErP0,PErPB]]] = state_info;
        let section = 'sub';
        for (var i = 0; i < M.length; i++) {
            if (area.x[i] <= 0) {
                section = 'sub';
            } else {
                section = 'sup';
            }
            M[i] = AMR(area.a[i]/At,g,section);
            P[i] = Pc*Math.pow(1+M[i]*M[i]*(g-1)/2,-g/(g-1));
            T[i] = Tc/(1+M[i]*M[i]*(g-1)/2);
            V[i] = M[i]*Math.sqrt(g*R*T[i]);
        }

    } else if (state_info[0] == 'subsonic') {
        let [MEsub, [PErP0sub,PErPB]] = state_info;
        for (var i = 0; i < M.length; i++) {
            M[i] = AMR(area.a[i]/At,g,'sub');
            P[i] = Pc*Math.pow(1+M[i]*M[i]*(g-1)/2,-g/(g-1));
            T[i] = Tc/(1+M[i]*M[i]*(g-1)/2);
            V[i] = M[i]*Math.sqrt(g*R*T[i]);
        }

    } else {
        let [state,[[M1,M2,ME], [P1rP01,P2rP1,PErP0,PErPB,P02rP01], [ANSErAt,AT2rAT]]] = state_info;
        let AreaShock = ANSErAt*At;
        let P02 = P02rP01*Pc;
        let At2 = AT2rAT*At;
        let section = 'sub';
        for (var i = 0; i < M.length; i++) {
            if ((area.x[i] <= 0)||(area.a[i]>AreaShock)) {
                section = 'sub';
            } else {
                section = 'sup';
            }
            if ((area.x[i] > 0)&&(area.a[i]>AreaShock)) {
                M[i] = AMR(area.a[i]/At2,g,section);
                P[i] = P02*Math.pow(1+M[i]*M[i]*(g-1)/2,-g/(g-1));
            } else {
                M[i] = AMR(area.a[i]/At,g,section);
                P[i] = Pc*Math.pow(1+M[i]*M[i]*(g-1)/2,-g/(g-1));
            }
            T[i] = Tc/(1+M[i]*M[i]*(g-1)/2);
            V[i] = M[i]*Math.sqrt(g*R*T[i]);
        }
    }

    let result = {x: area.x, r: area.r, m: M, p: P, t: T, v: V};
    return result;
}

export function get_results(g,E,At,P0,Patm,type,Lf,Tc,R,mode) {
    if (!Lf) {
        var Lf = 1;
    }
    let state = find_state(P0,Patm,g,At,At*E);
    let vectors = write_vectors(state,g,At,E,Lf,type,P0,Tc,R);
    let nozzle_xy = format_data(vectors.x,vectors.r);
    let nozzle_m = format_data(vectors.x,vectors.m);
    let nozzle_p = format_data(vectors.x,vectors.p,1,0.001);
    let nozzle_t = format_data(vectors.x,vectors.t);
    let nozzle_v = format_data(vectors.x,vectors.v);
    let angles = [0.261799387,0.261799387]
    if (type == 'bell') {
        angles = bell_angle(E,Lf);
    }
    let result;
    if (mode == 1) {
        result = [nozzle_xy, nozzle_m, nozzle_p, nozzle_t, nozzle_v];
    } else {
        result = [nozzle_xy, state, nozzle_m,nozzle_p, nozzle_t, nozzle_v,angles];
    }
    
    return result;
}

// main script

let g = 1.4, E = 10, At = 0.05, P0 = 5000000, Patm = 100000, R = 300, Tc = 3000, type = 'conical', Lf = 1;

let [nozzle_xy, nozzle_m, nozzle_p, nozzle_t, nozzle_v] = get_results(g,E,At,P0,Patm,type,Lf,Tc,R,1);

export {nozzle_xy, nozzle_m, nozzle_p, nozzle_t, nozzle_v}