use framework::vector::{size2,uint2,int2};

pub type Mesh = size2;
pub const fn N(M:Mesh) -> algebra::Len { (M.x*M.y) as algebra::Len }
fn div_remu(n : u32, d : u32) -> (u32, u32) { (n/d, n%d) }
fn div_rem(n : i32, d : u32) -> (i32, i32) { (n/d as i32, n%d as i32) }
use algebra::Idx;
pub fn mesh<const M:Mesh>(i:Idx) -> uint2 { div_remu(i as u32,M.x).into() }
pub fn dmesh<const M:Mesh>(i:Idx,j:Idx) -> int2 { div_rem(j as i32 - i as i32, M.x).into() }

use crate::algebra;
pub type Field<T=f32,const M:Mesh> = algebra::Array<T,{N(M)}>;
pub type Matrix<'a, const M:Mesh> = crate::compose::BoxFn<'a,(Idx, Idx),f32>;
pub type Operator<T,const M:Mesh> = Box<dyn algebra::Operator<T,{N(M)}>>;
pub type OperatorOnce<T,const M:Mesh> = Box<dyn algebra::OperatorOnce<T,{N(M)}>>;
pub struct Equation<T=f32,const M:Mesh> {pub A : algebra::LU<{N(M)}>, pub B : Box<dyn Fn()->Box<OperatorOnce<T,M>>>} // Ax = Bx' + ...

fn matrix<F:Fn(uint2, int2)->f32, const M:Mesh>(f : F) -> impl Fn(Idx,Idx)->f32 { move |i,j| { f(mesh::<M>(i), dmesh::<M>(i,j)) } }

use framework::{core::{mask,sign,abs,sq}, vector::{xy,vec2}};
pub type BC = fn(u32,u32,i32)->f32;

pub trait Operators<const M:Mesh> {
    /*fn operator<T:std::iter::Sum,const f:fn(uint2,int2)->f32>(v:&dyn algebra::Vector<T,{M(N)}>) -> algebra::BoxVector<T,M> where f32:std::ops::Mul<T,Output=T> {
        algebra::BoxVector::new(algebra::matrix_mul::<_,_,_,{N(M)}>(matrix::<_,M>(&f),v))
    }
    fn eq<T:std::iter::Sum,A:Fn(uint2,int2)->f32, const B:fn(uint2,int2)->f32>(A:A) -> Equation<T,M> where f32:std::ops::Mul<T,Output=T> {
        Equation{ A: algebra::LU::new(matrix::<_,M>(A)), B: Self::operator::<B> }
    }*/
    fn operator<T:std::iter::Sum,F:Fn(uint2,int2)->f32+Clone+'static>(f : F) -> Operator<T,M> where f32:std::ops::Mul<T,Output=T> {
        box move |v| algebra::BoxVector::new(algebra::matrix_mul::<_,_,_,{N(M)}>(matrix::<_,M>(f.clone()),v))
    }
    fn operator_once<T:std::iter::Sum,F:Fn(uint2,int2)->f32+'static>(f : F) -> OperatorOnce<T,M> where f32:std::ops::Mul<T,Output=T> {
        box move |v| algebra::BoxVector::new(algebra::matrix_mul::<_,_,_,{N(M)}>(matrix::<_,M>(f),v))
    }
    fn eq<T:std::iter::Sum,A:Fn(uint2,int2)->f32, B:Fn()->Box<dyn Fn(uint2,int2)->f32>+'static>(A:A, B:B) -> Equation<T,M> where f32:std::ops::Mul<T,Output=T> {
        Equation{ A: algebra::LU::new(matrix::<_,M>(A)), B: box move ||box Self::operator_once(B()) }
    }

    //fn I(_:uint2, d:int2) -> f32 { mask(d==0, 1.) }
    fn I()->crate::compose::BoxFn<'static,(uint2,int2),f32>{crate::compose::BoxFn::new(|_,d| mask(d==0, 1.) )}
    fn border(xy{x,y}:uint2) -> bool { x==0 || x==M.x-1 || y==0 || y==M.y-1 }
    fn interior(p:uint2,predicate:bool,value:f32) -> f32 { mask(predicate && !Self::border(p), value) }
    //fn P(p:uint2, d:int2) -> f32 { Self::interior(p, d==0, 1.) }
    fn P()->crate::compose::BoxFn<'static,(uint2,int2),f32>{crate::compose::BoxFn::new(|p,d| Self::interior(p, d==0, 1.) )}

    #[allow(non_upper_case_globals)] const δ : vec2 = framework::vector::div_f32(1f32, M.as_f32());
    fn Dx(p:uint2, d:int2) -> f32 { Self::interior(p, (abs(d.x),d.y) == (1,0), sign(d.x) as f32/(2.*Self::δ.x)) } // ∂x
    fn Dy(p:uint2, d:int2) -> f32 { Self::interior(p, (d.x,abs(d.y)) == (0,1), sign(d.y) as f32/(2.*Self::δ.y)) } // ∂y
    //fn Δ(p:uint2, d:int2) -> f32 {
    fn Δ()->crate::compose::BoxFn<'static,(uint2,int2),f32>{crate::compose::BoxFn::new(|p,d:int2|{
        Self::interior(p, true, {
            match (abs(d.x),abs(d.y)) {
                (0,0) => -2.*(1./sq(Self::δ.x)+1./sq(Self::δ.y)),
                (1,0) => 1./sq(Self::δ.x),
                (0,1) => 1./sq(Self::δ.y),
                _ => 0.
            }
        })
    })}

    fn BC<const X: BC, const Y: BC>(p:uint2, d:int2) -> f32 {
        if p.x==0 || p.x==M.x-1 { X(M.x,p.x,d.x) } // Horizontal boundary condition kernel (and corners)
        else if p.y==0 || p.y==M.y-1 { Y(M.y,p.y,d.y) }
        else { 0. }
    }
}
