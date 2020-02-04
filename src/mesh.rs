#![allow(non_snake_case)]
use framework::vector::{int2,uint2,size2}; use crate::algebra::{self,Idx};
pub type Mesh = size2;
pub const fn N(M:Mesh) -> algebra::Len { (M.x*M.y) as algebra::Len }
fn div_remu(n : u32, d : u32) -> (u32, u32) { (n/d, n%d) }
fn div_rem(n : i32, d : u32) -> (i32, i32) { (n/d as i32, n%d as i32) }
pub fn mesh<const M:Mesh>(i:Idx) -> uint2 { div_remu(i as u32,M.x).into() }
pub fn dmesh<const M:Mesh>(i:Idx,j:Idx) -> int2 { div_rem(j as i32 - i as i32, M.x).into() }

pub type Field<T=f32,const M:Mesh> = algebra::Array<T,{N(M)}>;
pub fn map<T,F:Fn(uint2)->T,const M:Mesh>(f: F) -> Field<T,M> { algebra::map(|i|f(mesh::<M>(i))) }

pub type Matrix<'a, const M:Mesh> = crate::compose::BoxFn<'a,(uint2,int2),f32>;
fn matrix<F:Fn(uint2,int2)->f32, const M:Mesh>(f : F) -> impl Fn(Idx,Idx)->f32 { move |i,j| { f(mesh::<M>(i), dmesh::<M>(i,j)) } }

pub type Operator<T,const M:Mesh> = Box<dyn algebra::Operator<T,{N(M)}>>;
pub fn operator<T:algebra::Sum,F:Fn(uint2,int2)->f32+Clone+'static,const M:Mesh>(f : F) -> Operator<T,M> where f32:std::ops::Mul<T,Output=T> {
    box move |v| algebra::BoxVector::new(algebra::matrix_mul::<_,_,_,{N(M)}>(matrix::<_,M>(f.clone()),v))
}

pub type BoxOperatorOnce<T,const M:Mesh> = Box<dyn algebra::OperatorOnce<T,{N(M)}>>;
fn box_operator_once<T:algebra::Sum, F:Fn(uint2,int2)->f32+'static, const M:Mesh>(f : F) -> BoxOperatorOnce<T,M> where f32:std::ops::Mul<T> { //{,Output=T> {
    box move |v| algebra::BoxVector::new(algebra::matrix_mul::<_,_,_,{N(M)}>(matrix::<_,M>(f),v))
}

//pub struct Equation<T=f32,const M:Mesh> {pub A : Box<algebra::LU<{N(M)}>>, pub B : Box<dyn Fn()->BoxOperatorOnce<T,M>>} // Ax = Bx' + ...
pub struct Equation<T=f32,const M:Mesh> {pub A : algebra::LU, pub B : Box<dyn Fn()->BoxOperatorOnce<T,M>>} // Ax = Bx' + ...
impl<T:algebra::Sum,const M:Mesh> Equation<T,M> where f32:std::ops::Mul<T> {
    pub fn new<A:Fn(uint2,int2)->f32, B:Fn()->Box<dyn Fn(uint2,int2)->f32>+'static>(A:A, B:B) -> Self {
        //Self{ A: algebra::LU::new(matrix::<_,M>(A)), B: box move ||box_operator_once::<_,_,M>(B()) }
        //Self{ A: algebra::LU::new::<_,{N(M)}>(matrix::<_,M>(A)), B: box move ||box_operator_once::<_,_,M>(B()) }
        Self{ A: algebra::LU::new(N(M), matrix::<_,M>(A)), B: box move ||box_operator_once::<_,_,M>(B()) }
    }
}

use framework::{core::{mask,sign,abs,sq}, vector::xy};
pub fn I<const M:Mesh>(_:uint2, d:int2) -> f32 { mask(d==0, 1.) }
pub fn border<const M:Mesh>(xy{x,y}:uint2) -> bool { x==0 || x==M.x-1 || y==0 || y==M.y-1 }
fn interior<const M:Mesh>(p:uint2,predicate:bool,value:f32) -> f32 { mask(predicate && !border::<M>(p), value) }
pub fn P<const M:Mesh>(p:uint2, d:int2) -> f32 { interior::<M>(p, d==0, 1.) }
fn δ<const M:Mesh>() -> framework::vector::vec2 { framework::vector::div_f32(1f32, M.as_f32()) }
pub fn Dx<const M:Mesh>(p:uint2, d:int2) -> f32 { interior::<M>(p, (abs(d.x),d.y) == (1,0), sign(d.x) as f32/(2.*δ::<M>().x)) } // ∂x
pub fn Dy<const M:Mesh>(p:uint2, d:int2) -> f32 { interior::<M>(p, (d.x,abs(d.y)) == (0,1), sign(d.y) as f32/(2.*δ::<M>().y)) } // ∂y
pub fn Δ<const M:Mesh>(p:uint2, d:int2) -> f32 {
    interior::<M>(p, true, {
    match (abs(d.x),abs(d.y)) {
            (0,0) => -2.*(1./sq(δ::<M>().x)+1./sq(δ::<M>().y)),
            (1,0) => 1./sq(δ::<M>().x),
            (0,1) => 1./sq(δ::<M>().y),
            _ => 0.
        }
    })
}
macro_rules! matrix_new_Op_M { ($($Op:ident)+) => ($( #[macro_export] macro_rules! $Op { () => ( Matrix::new($Op::<M>) ) } )+) } matrix_new_Op_M!(I P Δ);
macro_rules! operator_Op { ($($Op:ident)+) => ($( #[macro_export] macro_rules! $Op { () => ( operator::<_,_,M>($Op::<M>) ) } )+) } operator_Op!(Dx Dy);

//pub fn BC<const M:Mesh, const X: BC, const Y: BC>(p:uint2, d:int2) -> f32 {
#[macro_export] macro_rules! BC { ($X:ident, $Y:ident) => (
    {fn BC<const M:Mesh>(p:uint2, d:int2) -> f32 {
        if p.x==0 || p.x==M.x-1 { $X::<{M.x}>(p.x,d.x) } // Horizontal boundary condition kernel (and corners)
        else if p.y==0 || p.y==M.y-1 { $Y::<{M.y}>(p.y,d.y) }
        else { 0. }
    } BC::<M> }
)}
#[macro_export] macro_rules! operator_BC { ($X:ident, $Y:ident) => ( operator::<_,_,M>(BC!($X,$Y)) ) }
