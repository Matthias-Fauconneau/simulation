#![allow(incomplete_features,non_snake_case)]#![feature(const_generics,non_ascii_idents,box_syntax)]use simulation::*;

struct Quantity<T=f32,const M:Mesh> { u : Field<T,M>, f : [Field<T,M>; 2] }
impl<T:Zero, const M:Mesh> Quantity<T,M> { fn new<U0:Fn(uint2)->T>(u0: U0) -> Self { Self{u:mesh::map(u0), f:Zero::zero() } } }
impl<T:Copy+Add+Sum, const M:Mesh> Quantity<T,M> where f32:Mul<T>, LU:Solve<T> {
    fn step<G:Fn(Idx)->T>(&mut self, equation:&Equation<M>, δt:f32, g:&G) where T:Sub {
        // 2-step Adams-Bashforth:  y[2] = y[1] + 3/2·δt·f(t[1],y[1]) - 1/2·δt·f(t[0],y[0])
        let time = std::time::Instant::now();
        fn mul(a:f32, b:f32) -> f32 { a*b }
        self.u = equation.A.solve(|i| dot(&(equation.B)(i), &self.u) + (mul(3.,δt)/2.)*self.f[1][i] - (δt/2.)*self.f[0][i] + g(i));
        log!("solve", time.elapsed().as_millis());
    }
    fn advect(&mut self, v : &xy<Field<f32,M>>) {
        self.f.swap(0, 1);
        self.f[1] = map(|i| v.x[i]*dot(&Dx!()(i), &self.u) + v.y[i]*dot(&Dy!()(i), &self.u));
     }
}

use boundary_condition::*;
mod boundary_condition {
    use framework::core::{Zero,array::map};
    pub fn constant(_:u32) -> [(i32, f32); 3] { [(0, 1.),Zero::zero(),Zero::zero()] }
    fn kernel<const K:[f32;3], const C:[f32;2]>(p:u32) -> [(i32, f32); 3] {
        if p==0 { map(|i| (i as i32, C[0]*K[i])) } else { map(|i| (-(i as i32), C[1]*K[i])) }
    }
    pub fn derivative(p:u32) -> [(i32, f32); 3] { kernel::<{[-3.,4.,-1f32]},{[1./2.,-1./2.]}>(p) } // Boundary condition on derivative
    pub fn Thom(p:u32) -> [(i32, f32); 3] { kernel::<{[0.,-8.,1f32]},{[1.,1.]}>(p) } // Thom boundary condition
}

struct System<const M:Mesh> {
    δt : f32,
    T : Equation<M>, // Temperature
    ω : Equation<M>, // Vorticity
        ωT : f32, // Boussinesq approximation in buoyancy-driven flows
    φ : Equation<M>, // Stream function (u=∇×φ)
    C : Equation<M>, // Color (visualization)
}
impl<const M:Mesh> System<M> {
    fn new(δt : f32, Pr : f32, Ra : f32) -> Self {
        Self{δt,
            T: Equation::new(P!() - (    δt/2.)*Δ!() + BC!(constant,derivative),     P!() + (     δt/2.)*Δ!()),
            ω: Equation::new(I!() - (Pr*δt/2.)*Δ!()                                          ,     P!() + (Pr*δt/2.)*Δ!()), ωT: Ra*Pr/2.,
            φ: Equation::new(I!()-P!()          + Δ!()                                         ,-1.*P!() + (         0.)*Δ!()),
            C: Equation::new(P!() - (Pr*δt/2.)*Δ!() + BC!(constant,constant),     P!() + (Pr*δt/2.)*Δ!()),
        }
    }
}

struct State<const M:Mesh> {
    φ: [Field<f32,M>; 2],
    T : Quantity<f32,M>,
    ω : Quantity<f32,M>,
    C : Quantity<v2<f32>,M>
}

impl<const M:Mesh> State<M> {
    fn C0(p:uint2) -> v2<f32> { let v = p.as_f32() / (M - 1.into()).as_f32(); v2(v.x, v.y) }
    fn new() -> Self { Self{φ:Zero::zero(), T:Quantity::new(|_|0.), ω:Quantity::new(|_|0.), C:Quantity::new(Self::C0) } }
    fn C_BC(i:Idx) -> v2<f32> { let p = mesh::<M>(i); mask(border::<M>(p), Self::C0(p)) }
    fn step(&mut self, system : &System<M>) {
        // Solves implicit evolution
        let Self{φ, T, ω, C} = self;
        let mut T0 : Field<f32,M> = T.u.clone();
        let time = std::time::Instant::now();
        T.step(&system.T, system.δt, &|_|0.);
        log!("T", time.elapsed().as_millis()); let time = std::time::Instant::now();
        for i in 0..mesh::N(M) { T0[i] = system.ωT*(T0[i]+T.u[i]); } let T12=T0;
        let mut φp : Field<f32,M> = Zero::zero(); // fixme: may be uninitialized
        for i in 0..mesh::N(M) { φp[i] = 2.*φ[1][i]-φ[0][i]; }
        ω.step(&system.ω, system.δt, &|i| dot(&Dx!()(i),&T12) + dot(&rows_BC!(Thom,Thom)(i), &φp));
        log!("w", time.elapsed().as_millis()); let time = std::time::Instant::now();
        C.step(&system.C, system.δt, &Self::C_BC);
        log!("C", time.elapsed().as_millis()); let time = std::time::Instant::now();
        // Evaluates explicit advection
        φ.swap(0, 1); φ[1]=system.φ.A.solve(|i|dot(&(system.φ.B)(i),&ω.u));
        log!("φ", time.elapsed().as_millis()); let time = std::time::Instant::now();
        let v = &xy{x:map(|i|dot(&Dy!()(i),&φ[1])), y: map(|i|-dot(&Dx!()(i),&φ[1]))};
        log!("v", time.elapsed().as_millis()); let time = std::time::Instant::now();
        T.advect(v);
        log!("T advect", time.elapsed().as_millis()); let time = std::time::Instant::now();
        ω.advect(v);
        log!("w advect", time.elapsed().as_millis()); let time = std::time::Instant::now();
        C.advect(v);
        log!("explicit C", time.elapsed().as_millis());
   }
}

struct Simulation<const M:Mesh> {
    system : System<M>,
    state : State<M>
}
impl<const M:Mesh> Simulation<M> { fn new(δt : f32, Pr : f32, Ra : f32) -> Self { Self{system: System::new(δt, Pr, Ra), state: State::new()} } }

impl<const M:Mesh> Solution for Simulation<M> {
    fn current(&self) -> view::Field { view::Field{size: M, values: &self.state.T.u.0} }
    fn step(&mut self) { self.state.step(&self.system); }
}

fn main() -> Result {
    const M:Mesh = xy{x:128,y:128};
    struct Parameters {Pr : f32, Ra : f32}
    fn parameters() -> Parameters {
        use SI::*;
        let η : DynamicViscosity = 2e-5 |Pa·s; //kg/m/s
        let ρ : MassDensity = 1.2 |kg_m3; //kg/m³
        let α : ThermalDiffusivity = 2e-5 |m2_s; //m²/s
        let ΔT : Temperature = 1.0 |K; // Difference between walls
        let β : ThermalExpansion = 1./(300.|K); //K¯¹: Ideal gas
        let g : Acceleration = 9.8 |m_s2; //m/s²: Gravity
        let L : Length = 0.1 |m; // Box length
        //let t : Time = sq(L)/α; //L²/α s: Box thermal diffusion time
        let ν : KinematicViscosity = η/ρ; //m²/s
        Parameters{Pr: (ν/α).into(), //Prandtl: kinematic/thermal diffusivity
                            Ra: ((ΔT*β*g*cb(L))/(ν*α)).into() //ΔT·β·g·L³/(ν·α) Rayleigh: Laminar/turbulent flow
        }
    }
    let Parameters{Pr,Ra} = parameters();
    let δt : f32 = 1./(max(M.x,M.y) as f32*sqrt(Ra)); //R·√Ra 1: Temporal resolution
    run(View!(Simulation::<M>::new(δt, Pr, Ra)))
}
