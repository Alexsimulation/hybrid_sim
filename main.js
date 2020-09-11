function ISP_F(of_in,props) {
    let pns = [0, 0, 0, 0, 0, 0, 0];
    let out_range = 0;
    if (props == 'N2O/Paraffin') {
        pns = [-0.0016, 0.0585, -0.8467, 6.0614, -23.0533, 58.0331, 179.7000];
        out_range = 300;
    } else if (props == 'LOX/Paraffin') {
        pns = [-0.0011, 0.0791, -1.8579, 20.4487, -114.4732, 296.5522, 69.1333];
        out_range = 255;
    } else if (props == 'LOX/PMMA') {
        pns = [-0.0107, 0.4045, -6.0904, 46.2984, -184.2097, 337.3317, 96.4000];
        out_range = 194;
    } else if (props == 'N2O/PMMA') {
        pns = [-0.0005, -0.0014, 0.3253, -4.2876, 16.8540, 3.5648, 198.5667];
        out_range = 260;
    }
    let ISP = 0;
    if (of_in < 10) {
        let i = pns.length;
        while (i >= 1) {
            ISP = ISP + pns[pns.length-i]*Math.pow(of_in,i-1);
            i = i - 1;
        }
    } else {
        ISP = out_range*Math.pow(2.718,-0.2*(of_in-10));
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

    let ind = 0;

    while ((r < re)&(ind <= 1000)) {
        Gox = m_dot_ox/(Math.PI*r*r);
        r_dot = a*Math.pow(Gox,n);
        r = r + r_dot*dt;
        m_dot_fuel = rho_fuel*L*Math.PI*(2*r*r_dot + dt*r_dot*r_dot);
        of = m_dot_ox/m_dot_fuel;
        T = T_F(m_dot_ox,m_dot_fuel,of,propellants);
        I = I + T*dt;
        t = t + dt;

        radius_data = radius_data.concat({x: t, y: r*100});
        thrust_data = thrust_data.concat({x: t, y: T/1000});
        impulse_data = impulse_data.concat({x: t, y: I/1000});
        of_data = of_data.concat({x: t, y: of});
        
        ind = ind + 1;
    }
    let values = [radius_data, thrust_data, impulse_data, of_data];
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

let [radius_data, thrust_data, impulse_data, of_data] = get_results(L,ri,re,m_dot_ox,propellants,a,n);

export {radius_data, thrust_data, impulse_data, of_data}