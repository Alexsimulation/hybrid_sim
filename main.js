function ISP_F(x) {
    let pns = [-0.0016, 0.0585, -0.8467, 6.0614, -23.0533, 58.0331, 179.7000];
    let ISP = 0;
    let i = pns.length;
    while (i >= 1) {
        ISP = ISP + pns[pns.length-i]*Math.pow(x,i-1);
        i = i - 1;
    }
    return ISP;
}
function T_F(m1,m2,x1){
    let T_1 = (m1+m2)*ISP_F(x1)*9.81;
    return T_1;
}
function get_rho(fuels) {
    let rho = 1000;
    if (fuels == 'N2O/Paraffin') {
        rho = 900;
    } else if (fuels == 'NaN') {
        rho = 1000;
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
    let T = T_F(m_dot_ox,m_dot_fuel,of);
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
        T = T_F(m_dot_ox,m_dot_fuel,of);
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

let ri = 0.2;
let re = 0.3;
let L = 1;
let m_dot_ox = 2;
let propellants = 'N2O/Paraffin';

let a = 0.00015;
let n = 0.5;

let [radius_data, thrust_data, impulse_data, of_data] = get_results(L,ri,re,m_dot_ox,propellants,a,n);

export {radius_data, thrust_data, impulse_data, of_data}