use framework::{core::{abs,mask}, vector::{size2,uint2,int2}};
use crate::algebra;

pub type Mesh = size2;
pub type Operator<const M:Mesh> = Box<dyn algebra::Matrix<{M.x*M.y}>>;
pub struct Equation<const M:Mesh> { A : algebra::LU<{M.x*M.y}>, B : Operator<M> } // Ax = Bx' + ...
type DenseField<T=f32,const M:Mesh> = algebra::DenseVector<T,{M.x*M.y}>;
pub type Field<T=f32,const M:Mesh> = DenseField<T,M>;

fn div_remu(n : u32, d : u32) -> (u32, u32) { (n/d, n%d) }
fn div_rem(n : i32, d : u32) -> (i32, i32) { (n/d as i32, n%d as i32) }
fn map<T,F:Fn(uint2)->T,const M:Mesh>(f : F) -> impl algebra::SizedVector<T,{M.x*M.y}> { move |i|{ f(div_remu(i,M.x).into()) } }
pub fn field<T,F:Fn(uint2)->T,const M:Mesh>(f : F) -> algebra::DenseVector<T,{M.x*M.y}> { algebra::collect(map(M,f)) }
//fn operator<F:Fn(uint2, int2)->f32,const N:u32>(M:size2, f : F) -> impl algebra::Matrix<N> { move |i,j| { f(div_remu(i,M.x).into(), div_rem(j as i32 - i as i32, M.x).into()) } }
//fn op<'a,F:Fn(uint2,int2)->f32+'a, const N:u32>(M:size2, f : F) -> framework::compose::RcFn<'a,(u32, u32),f32> { framework::compose::RcFn::new::<F>(operator(M,f)) }
pub fn op<F:Fn(uint2,int2)->f32+'static, const M:Mesh>(f : F) -> framework::compose::RcFn<'static,(u32, u32),f32> { 
    framework::compose::RcFn::new( move |i,j| { f(div_remu(i,M.x).into(), div_rem(j as i32 - i as i32, M.x).into()) } )
}
pub fn eq<A:algebra::Matrix<M>, B:algebra::Matrix<M>+'static, const M:Mesh>(A:A, B:B) -> Equation<M> { Equation{ A: algebra::LU::new(A), B: box B } }

pub fn identity(d:int2) -> f32 { mask(d==0, 1.) }
fn BC(sym:f32,c:[f32;3],M:u32,p:u32,d:i32) -> f32 { // Boundary condition kernel
    if (-(c.len() as i32)+1..c.len() as i32).contains(&d) {
        if p==0 { return c[d as usize] }
        else if p==M-1 { return sym*c[(M-1-abs(d) as u32) as usize] }
    }
    0.
}
pub fn BCx(sym:f32,c:[f32;3],M:size2,p:uint2,d:int2) -> f32 {                                                    BC(sym,c,M.x,p.x,d.x)/(1. / M.x as f32) } // Horizontal boundary condition kernel
pub fn BCy(sym:f32,c:[f32;3],M:size2,p:uint2,d:int2) -> f32 { mask((1..M.x-1).contains(&p.x), BC(sym,c,M.y,p.y,d.y)/(1. / M.y as f32))} // Corners use BCx
