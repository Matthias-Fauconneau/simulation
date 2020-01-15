#![allow(incomplete_features,uncommon_codepoints)]#![feature(const_generics,non_ascii_idents,fn_traits,unboxed_closures,trait_alias,box_syntax)]
//use std::mem::swap; 
use framework::{core::{Zero,mask,sign,abs,sq,cb,sqrt}, vector::{xy,uint2,vec2}};
mod algebra;
mod mesh; use mesh::{Mesh, Operator,Equation,Field, op,eq,field, identity,BCx,BCy};

struct System<const M:Mesh> {
    D: xy<Operator<M>>,
    T : Equation<M>, // Temperature
    ω : Equation<M>, ωT : Operator<M>, ωφ : Operator<M>, // Vorticity
    φ : Equation<M>, // Stream function (u=∇×φ)
    𝜓 : Equation<M>, 𝜓_G : Box<dyn Fn(u32)->(f32,f32)>, // Color (visualization)
}

impl<const M:Mesh> System<M> {
    fn new(δt : f32, Pr : f32, Ra : f32) -> Self {
        let I = &op::<_,M>(|_,d| { identity(d) });
        let border = |M,xy{x,y}| -> bool { x==0 || x==M.x-1 || y==0 || y==M.y-1 };
        let interior = |M,p:uint2,predicate:bool,value:f32| -> f32 { mask(predicate && !border(M,p), value) };
        let P = &op::<_,M>(|p,d| { interior(M,p, d==0, 1.) });
        
        let δ = vec2{x: 1./(M.x as f32), y: 1./(M.y as f32) };  // : vec2 = 1f32/M.into();
        let D = xy{ x:&op(M,|p,d| { interior(M,p, (abs(d.x),d.y) == (1,0), sign(d.x) as f32/(2.*δ.x)) }), // ∂x
                        y:&op(M,move |p,d| { interior(M,p, (d.x,abs(d.y)) == (0,1), sign(d.y) as f32/(2.*δ.y)) })}; // ∂y
        let Δ = &op(M,move |p,d| { interior(M,p, true, {
            match (abs(d.x),abs(d.y)) {
                (0,0) => -2.*(1./sq(δ.x)+1./sq(δ.y)),
                (1,0) => 1./sq(δ.x),
                (0,1) => 1./sq(δ.y),
                _ => 0.
            }
        })});
        
        let BC_T = op::<_,M>(|p,d| { if p.x==0 || p.x==M.x-1 { identity(d) } else { BCy(-1., [-3.,4.,-1.],M,p,d)/2. } }); // constant value on vertical, derivative on horizontal
        let ωφ = box op::<_,M>(|p,d| { let thom=[0.,-8.,1.]; (if p.x==0 || p.x==M.x-1 { BCx } else { BCy })(1.,thom,M,p,d) }); // Thom horizontal
        let BC_𝜓 = op::<_,M>(|p,d| { mask(border(M,p) && d==0, 1.) }); // Constant boundary condition for advection (coefficients)
        let 𝜓_G = box move |p| { mask(border(M,p), (p.x as f32/(M.x-1) as f32, p.y as f32/(M.y-1) as f32)) }; // Constant boundary condition for advection (= source constants)
        
        Self{D:xy{x:box D.x, y:box D.y}, //box D,
            T: eq(P      - (δt/2.)*Δ + BC_T  ,      P + δt/2.*Δ     ),
            ω: eq(I  - (Pr*δt/2.)*Δ               ,      P + Pr*δt/2.*Δ), ωT: box( Ra*Pr*δt/2.*D.x ), ωφ,
            φ:  eq(I-P             + Δ               ,-1.*P                    ),
            𝜓: eq(P - (Pr*δt/2.)*Δ + BC_𝜓,     P + Pr*δt/2.*Δ), 𝜓_G }
    }
}

struct Tω𝜓<const M:Mesh> { T : Field<f32,M>, ω : Field<f32,M>, 𝜓 : Field<vec2,M> } 
// ICE: Field<N>. const argument index assumes optional type argument was given and gets OOB ~ 67858
impl<const M:Mesh> Zero for Tω𝜓<M>{ fn zero() -> Self {Self{T:Zero::zero(),ω:Zero::zero(),𝜓:Zero::zero()}} } // fixme: #[derive(Zero)]
//impl<const N: u32> Zero for Tω𝜓<N>{ fn zero() -> Self {Self{T:Field::<f32,N>::zero(),ω:Field::<f32,N>::zero(),𝜓:Field::<vec2,N>::zero()}} } // fixme: Zero::zero()
//fn mul<T:std::ops::Mul<Field<f32,N>>, const N:u32>(a: T, b: &Tω𝜓<N>) -> Tω𝜓<N> where <T as std::ops::Mul<Field<f32,N>>>::Output:Into<Field<f32,N>> { 
//fn mul<const N:u32>(a: f32, b: &Tω𝜓<N>) -> Tω𝜓<N> { Tω𝜓::<N>{T: a*b.T, ω: a*b.ω, 𝜓: a*b.𝜓} }
//fn mul<const N:u32>(a: f32, b: &Tω𝜓<N>) -> Tω𝜓<N> { Tω𝜓::<N>{T: a*b.T, ω: a*b.ω, 𝜓: a*operator::RcFn::new(&b.𝜓)} }
//impl<const N:u32> std::ops::Mul<&Tω𝜓<N>> for Op<'_> { type Output=Tω𝜓<N>; fn mul(self, b: &Tω𝜓<N>) -> Self::Output { mul(self, b) } }
//impl<const N:u32> std::ops::Mul<&mut Tω𝜓<N>> for Op<'_> { type Output=Tω𝜓<N>; fn mul(self, b: &mut Tω𝜓<N>) -> Self::Output { mul(self, b) } }

#[allow(non_camel_case_types)] struct φA<const M:Mesh> { φ: Field<M>, A: Tω𝜓<M> }
impl<const M:Mesh> Zero for φA<M> { fn zero() -> Self {Self{φ:Field::<f32,M>::zero(), A:Zero::zero()}} } // fixme: #[derive(Zero)]
struct State<const M:Mesh> {
    φA : [φA<M>; 2], // [previous,current] {stream function, non-linear advection term}
    C : Tω𝜓<M>,
}
impl<const M:Mesh> Zero for State<M> { fn zero() -> Self {Self{φA:Zero::zero(), C:Zero::zero()}} } // fixme: #[derive(Zero)]
impl<const M:Mesh> State<M> {
    fn new() -> Self { Self{C:Tω𝜓{ 𝜓:field(M, |p:uint2|->vec2 { p.as_f32() / (M-1.into()).as_f32() }), ..Zero::zero() }, ..Zero::zero()} }
    fn update(&mut self, system : &System<M>, _δt : f32) {
        let _A = system.T.B*self.C.T;
        //let _A = system.T.B(self.C.T);
        //framework::core::log( A );
        // Solves implicit evolution
        /*let &System{D, T, ω, ωT, ωφ, φ, 𝜓, 𝜓_G} = &system;
        let &mut Self{C, φA:[φA{φ:φp, A:Ap}, φA{φ:φc, A:Ac}]} = &mut self;
        use std::ops::Mul; let B : operator::BoxFn<(u32,),f32> = 3f32.mul(Ac.T);
        let B = 3f32*Ac.T;
        //let B : operator::BoxFn<'_,(u32,),f32> = 3f32*Ac.T;
        //let test = A + B;
        let test = T.B*C.T + 3f32*Ac.T;
        let Tn = T.A.solve(T.B*C.T  + (3.*δt/2.)*Ac.T - (δt/2.)*Ap.T);
        *C=Tω𝜓{ ω: ω.A.solve(ω.B*C.ω + (3.*δt/2.)*Ac.ω - (δt/2.)*Ap.ω + ωT(C.T+&Tn) + ωφ*(2.*φc-φp)),
                        𝜓: 𝜓.A.solve(𝜓.B*C.𝜓 + (3.*δt/2.)*Ap.𝜓  - (δt/2.)*Ap.𝜓  + 𝜓_G), T: Tn};
        // Evaluates explicit advection
        swap(φp, φc); *φc=φ.A.solve(φ.B*C.ω);
        let U = xy{x:D.y, y:-D.x}*φc; //cross(D)*φc; // need cross = complicated
        swap(Ap, Ac); *Ac = (U.x*D.x + U.y*D.y)*C; // dot(U,  D)*C; // need dot by ref (+lifetime) = complicated*/
    }
}

struct Parameters {Pr : f32, Ra : f32}
mod SI;
fn parameters() -> Parameters {
    use crate::SI::*;
    let η : DynamicViscosity = 2e-5 |Pa·s; //kg/m/s
    let ρ : MassDensity = 1.2 |kg_m3; //kg/m³
    let α : ThermalDiffusivity = 2e-5 |m2_s; //m²/s
    let ΔT : Temperature = 1.0 |K; // Difference between walls
    let β : ThermalExpansion = 1./(300.|K); //K¯¹: Ideal gas
    let g : Acceleration = 9.8 |m_s2; //m/s²: Gravity
    let L : Length = 0.1 |m; // Box side
    //let t : Time = sq(L)/α; //L²/α s: Box thermal diffusion time
    let ν : KinematicViscosity = η/ρ; //m²/s
    Parameters{Pr: ν/α, //Prandtl: kinematic/thermal diffusivity
                        Ra: (ΔT*β*g*cb(L))/(ν*α) //ΔT·β·g·L³/(ν·α) Rayleigh: Laminar/turbulent flow
    }
}

fn main() {
    pub const R : u32 = 128; // Resolution
    pub const M : uint2 = xy{x:R+1, y:R+1};
    let Parameters{Pr,Ra} = parameters();
    let δt : f32 = 1./(R as f32*sqrt(Ra)); //R·√Ra 1: Temporal resolution
    let system = System::<M>::new(δt, Pr, Ra);
    let mut state = State::new();
    state.update(&system, δt);
    /* subplot(position, size, 4, 0, Cw, Mx, My, "Vorticity ω"_);
        subplot(position, size, 4, 1, Ux, Uy, Mx, My, "Velocity u"_);
        subplot(position, size, 4, 2, positiveToImage(Ct, Mx, My), "Temperature T"_);
        subplot(position, size, 4, 3, toImage(CTx, CTy, Mx, My), "Advection"_);*/
}
