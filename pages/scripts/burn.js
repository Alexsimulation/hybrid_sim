function ISP_F(of_in,props) {
    let pns = [0, 0, 0, 0, 0, 0, 0];
    let out_range = 0;
    let out_coef = 0;
    if (props == 'N2O/Paraffin') {
        pns = [-0.0015972222222,   0.0584775641026,  -0.8466613247863,   6.0614291958040,  -23.0533148795640,   58.0330652680635,   179.7000000000013];
        out_range = 300.04930069930122726873;
        out_coef = 0.01770377605488909073 ;
    } else if (props == 'LOX/Paraffin') {
        pns = [-0.0011111111111,   0.0791025641025,  -1.8579059829055,   20.4486596736565,  -114.4731507381392,   296.5521911421711,   69.1333333333463];
        out_range = 256.08531468531361952046;
        out_coef = 0.00903569618427984417;
    } else if (props == 'LOX/PMMA') {
        pns = [-0.0106944444444,   0.4044551282051,  -6.0904380341877,   46.2984120046597,  -184.2097066821982,   337.3317482517328,   96.4000000000109];
        out_range = 193.84685314684838886023;
        out_coef = 0.06265436119621256572;
    } else if (props == 'N2O/PMMA') {
        pns = [-0.0004861111111,  -0.0013942307692,   0.3253472222221,  -4.2875801282041,   16.8539850427316,   3.5647902097956,   198.5666666666635];
        out_range = 259.97097902097897303975;
        out_coef = 0.02037496999093071914;
    }
    let ISP = 0;
    if (of_in < 10) {
        let i = pns.length;
        while (i >= 1) {
            ISP = ISP + pns[pns.length-i]*Math.pow(of_in,i-1);
            i = i - 1;
        }
    } else {
        ISP = out_range*Math.pow(2.718,-out_coef*(of_in-10));
    }
    return ISP;
}
function T_F(m1,m2,x1,props){
    let T_1 = (m1+m2)*ISP_F(x1,props)*9.81;
    return T_1;
}
function get_rho(fuels) {
    let rho = 900;
    if (fuels == 'N2O/Paraffin') {
        rho = 900;
    } else if (fuels == 'LOX/Paraffin') {
        rho = 900;
    } else if (fuels == 'LOX/PMMA') {
        rho = 1170;
    } else if (fuels == 'N2O/PMMA') {
        rho = 1170;
    }
    return rho;
}
export function get_results(L,ri,re,m_dot_ox,propellants,a,n) {
    let rho_fuel = get_rho(propellants);
    let t = 0;
    let dt = (re-ri)/(a*Math.pow((4*m_dot_ox/(Math.PI*Math.pow((re+ri),2))),n)*200);
    let Gox = m_dot_ox/(Math.PI*ri*ri);
    let r = ri;
    let r_dot = a*Math.pow(Gox,n);
    let m_dot_fuel = rho_fuel*L*Math.PI*(2*r*r_dot + dt*dt*r_dot*r_dot);
    let of = m_dot_ox/m_dot_fuel;
    let T = T_F(m_dot_ox,m_dot_fuel,of,propellants);
    let I = 0;
    
    let radius_data = [];
    let thrust_data = [];
    let impulse_data = [];
    let of_data = [];
    let ISP_data = [];

    let ind = 0;
    let last_r_dot = r_dot;
    let last_T = T;

    while ((r < re)&(ind <= 1000)) {
        Gox = m_dot_ox/(Math.PI*r*r);
        r_dot = a*Math.pow(Gox,n);
        m_dot_fuel = rho_fuel*L*Math.PI*(2*r*r_dot + dt*r_dot*r_dot);
        of = m_dot_ox/m_dot_fuel;
        T = T_F(m_dot_ox,m_dot_fuel,of,propellants);
        r = r + (r_dot+last_r_dot)*dt/2;
        I = I + (T+last_T)*dt/2;
        t = t + dt;

        last_r_dot = r_dot;
        last_T = T;

        radius_data = radius_data.concat({x: t, y: r*100});
        thrust_data = thrust_data.concat({x: t, y: T/1000});
        impulse_data = impulse_data.concat({x: t, y: I/1000});
        of_data = of_data.concat({x: t, y: of});
        ISP_data = ISP_data.concat({x: t, y: ISP_F(of,propellants)});
        
        ind = ind + 1;
    }
    let values = [radius_data, thrust_data, impulse_data, of_data, ISP_data];
    return values;
}

// main script

let ri = 0.05;
let re = 0.1;
let L = 0.5;
let m_dot_ox = 2;
let propellants = 'N2O/Paraffin';

let a = 0.00015;
let n = 0.5;

let [radius_data, thrust_data, impulse_data, of_data, ISP_data] = get_results(L,ri,re,m_dot_ox,propellants,a,n);

export {radius_data, thrust_data, impulse_data, of_data, ISP_data}