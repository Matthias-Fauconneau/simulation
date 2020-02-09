#![allow(non_snake_case)]
use {framework::{core::{Zero,div_rem,array::Iterator},vector::{int2,uint2,size2}},crate::algebra::{self,Idx,Array}};
pub type Mesh = size2;
pub const fn N(M:Mesh) -> usize { (M.x*M.y) as usize }
pub fn mesh<const M:Mesh>(i:Idx) -> uint2 { let (div, rem) = div_rem(i as u32,M.x); xy{x:rem, y:div} }
fn dmesh<const M:Mesh>(p: uint2, d: int2) -> u32 { ((p.y as i32+d.y) as u32)*M.x+(p.x as i32+d.x) as u32 }

pub type Field<T=f32,const M:Mesh> = Array<T,{N(M)}>;
pub fn map<T,F:Fn(uint2)->T,const M:Mesh>(f: F) -> Field<T,M> { algebra::map(|i|f(mesh::<M>(i))) }

/*pub trait Op<const KA:usize, const KB:usize> { type Output; fn op<const OP:fn(f32,f32)->f32>(A : Row<KA>, B: Row<KB>) -> Self::Output; }
impl<const KA:usize, const KB:usize> Op<KA,KB> for Row<KA> {
    type Output = Row<{KA+KB}>;
    fn op<const OP:fn(f32,f32)->f32>(A : Row<KA>, B: Row<KB>) -> Self::Output {
        let mut R : Self::Output = Zero::zero(); // fixme: may be uninitialized
        let (mut a, mut b, mut r) = (0,0,0);
        loop {
            if a < A.len() {
                if b < B.len() {
                    if A[a].0 < B[b].0 { R[r] = A[a]; a+=1; r+=1; }
                    else if B[b].0 < A[a].0 { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                    else { R[r] = (A[a].0, OP(A[a].1,B[b].1)); a+=1; b+=1; r+=1; }
                } else {
                    while a < A.len() { R[r] = A[a]; a+=1; r+=1; }
                    break;
                }
            } else {
                while b < B.len() { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                break;
            }
        }
        if r < KA+KB/*R.len()*/ { R[r] = Zero::zero(); } // Zero terminate
        R
    }
}

pub type Row<const K:usize> = Array<(int2,f32), K>;
impl<const A:usize, const B:usize> std::ops::Add<Row<B>> for Row<A> where Self:Op<A,B> {
    type Output=<Self as Op<A,B>>::Output;
    fn add(self, b: Row<B>) -> Self::Output { fn add(a:f32,b:f32)->f32{a+b} Row::op::<add>(self, b) }
}
impl<const A:usize, const B:usize> std::ops::Sub<Row<B>> for Row<A> {
    type Output=<Self as Op<A,B>>::Output;
    fn sub(self, b: Row<B>) -> Self::Output where Self::Output: { fn sub(a:f32,b:f32)->f32{a-b} Row::op::<sub>(self,b) }
}*/

pub type Row<const K:usize> = Array<(int2,f32), K>;
impl<const KA:usize, const KB:usize> std::ops::Add<Row<KB>> for Row<KA> {
    type Output = Row<{KA+KB}>;
    fn add(self, B: Row<KB>) -> Self::Output {
        fn OP(a:f32,b:f32)->f32{a+b}
        let A = self;
        let mut R : Self::Output = Zero::zero(); // fixme: may be uninitialized
        let (mut a, mut b, mut r) = (0,0,0);
        loop {
            if a < A.len() {
                if b < B.len() {
                    if A[a].0 < B[b].0 { R[r] = A[a]; a+=1; r+=1; }
                    else if B[b].0 < A[a].0 { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                    else { R[r] = (A[a].0, OP(A[a].1,B[b].1)); a+=1; b+=1; r+=1; }
                } else {
                    while a < A.len() { R[r] = A[a]; a+=1; r+=1; }
                    break;
                }
            } else {
                while b < B.len() { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                break;
            }
        }
        if r < KA+KB/*R.len()*/ { R[r] = Zero::zero(); } // Zero terminate
        R
    }
}
impl<const KA:usize, const KB:usize> std::ops::Sub<Row<KB>> for Row<KA> {
    type Output = Row<{KA+KB}>;
    fn sub(self, B: Row<KB>) -> Self::Output where Self::Output: {
        fn OP(a:f32,b:f32)->f32{a-b}
        let A = self;
        let mut R : Self::Output = Zero::zero(); // fixme: may be uninitialized
        let (mut a, mut b, mut r) = (0,0,0);
        loop {
            if a < A.len() {
                if b < B.len() {
                    if A[a].0 < B[b].0 { R[r] = A[a]; a+=1; r+=1; }
                    else if B[b].0 < A[a].0 { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                    else { R[r] = (A[a].0, OP(A[a].1,B[b].1)); a+=1; b+=1; r+=1; }
                } else {
                    while a < A.len() { R[r] = A[a]; a+=1; r+=1; }
                    break;
                }
            } else {
                while b < B.len() { R[r] = (B[b].0, OP(0.,B[b].1)); b+=1; r+=1; }
                break;
            }
        }
        if r < KA+KB/*R.len()*/ { R[r] = Zero::zero(); } // Zero terminate
        R
    }
}


impl<const K:usize> std::ops::Mul<Row<K>> for f32 {
    type Output=Row<K>;
    fn mul(self, b: Row<K>) -> Self::Output { Array(Iterator::collect(b.iter().map(|(d,v)|(*d,self**v)))) }
}

pub type Rows<'t, const M:Mesh, const K:usize> = crate::compose::BoxFn<'t,(uint2,),Row<K>>;
pub fn rows<F:Fn(uint2)->Row<K>, const M:Mesh, const K:usize>(f : F) -> impl /*algebra::Rows<{N(M)},{K}>*/ Fn(Idx)->algebra::Row<K> {
    move |i| {
        let p = mesh::<M>(i);
        //Array(Iterator::collect( f(p).iter().map(|(d,v)| (dmesh::<M>(p,*d), *v)) ))
        algebra::map({let f=f(p); move |j| (dmesh::<M>(p,f[j].0), f[j].1)})
    }
}

pub struct Equation<const M:Mesh> {pub A : algebra::LU, pub B : Box<dyn algebra::Rows<{N(M)},6>>} // Ax = Bx' + ...
impl<const M:Mesh> Equation<M> {
    pub fn new<A:Fn(uint2)->Row<KA>, B:Fn(uint2)->Row<6>+'static, const KA:usize>(A:A, B:B) -> Self {
        Self{A: algebra::LU::new(/*N(M),*/ &algebra::columns::<_,KA,{N(M)},6>(rows::<_,M,KA>(A))), B: box rows::<_,M,6>(B)}
    }
}

use framework::{core::{mask,sq}, vector::xy};
pub fn I<const M:Mesh>(_:uint2) -> Row<1> { Array([(xy{x:0, y:0}, 1.)]) }
pub fn border<const M:Mesh>(xy{x,y}:uint2) -> bool { x==0 || x==M.x-1 || y==0 || y==M.y-1 }
fn interior<const M:Mesh, const K:usize>(p:uint2, kernel: Row<K>) -> Row<K> { mask(!border::<M>(p), kernel) }
pub fn P<const M:Mesh>(p:uint2) -> Row<1> { interior::<M,1>(p, Array([(xy{x:0, y:0}, 1.)])) }
fn δ<const M:Mesh>() -> framework::vector::vec2 { framework::vector::div_f32(1f32, M.as_f32()) }
pub fn Dx<const M:Mesh>(p:uint2) -> Row<2> { let c = 1./(2.*δ::<M>().x); interior::<M,2>(p, Array([(xy{x:-1, y:0}, -c), (xy{x:1, y:0}, c)])) } // ∂x
pub fn Dy<const M:Mesh>(p:uint2) -> Row<2> { let c = 1./(2.*δ::<M>().y); interior::<M,2>(p, Array([(xy{x:0, y:-1}, -c), (xy{x:0, y:1}, c)])) } // ∂y
pub fn Δ<const M:Mesh>(p:uint2) -> Row<5> {
    let c = 1./sq(δ::<M>());
    interior::<M,5>(p, Array([
                                  (xy{x:0, y:-1}, c.y),
    (xy{x:-1, y:0}, c.x), (xy{x:0, y:0}, -2.*(c.x+c.y)), (xy{x:1, y:0}, c.x),
                                  (xy{x:0, y:1}, c.y) ]))
}
macro_rules! Rows_new_Op_M { ($($Op:ident)+) => ($( #[macro_export] macro_rules! $Op { () => ( $crate::mesh::Rows::new($Op::<M>) ) } )+) } Rows_new_Op_M!(I P Δ);
macro_rules! rows_Op { ($($Op:ident)+) => ($( #[macro_export] macro_rules! $Op { () => ( $crate::mesh::rows::<_,M,2>($Op::<M>) ) } )+) } rows_Op!(Dx Dy);

/*type BC = fn(u32)->[(i32, f32); 3];
pub fn BC<const X: BC, const Y: BC, const M:Mesh>(p:uint2) -> Row<3> {
    if p.x==0 || p.x==M.x-1 { let t=X(p.x); algebra::map(|i| (xy{x:t[i].0,y:0}, t[i].1)) } // Horizontal boundary condition kernel (and corners)
    else if p.y==0 || p.y==M.y-1 { let t=Y(p.y); algebra::map(|i| (xy{x:0,y:t[i].0}, t[i].1)) }
    else { Zero::zero() }
}
#[macro_export] macro_rules! BC { ($X:ident, $Y:ident) => ( $crate::mesh::Rows::new(($crate::mesh::BC::<$X,$Y,M>)) ) }
#[macro_export] macro_rules! rows_BC { ($X:ident, $Y:ident) => ( $crate::mesh::rows::<_,M,3>($crate::mesh::BC::<$X,$Y,M>) ) }*/
#[macro_export] macro_rules! BCxy { ($X:ident, $Y:ident) => (
    {fn BCxy<const M:Mesh>(p:uint2) -> $crate::mesh::Row<3> {
        if p.x==0 || p.x==M.x-1 { let t=$X::<{1./(M.y as f32)},{1./(M.x as f32)}>(p.x); algebra::map(|i| (xy{x:t[i].0,y:0}, t[i].1)) } // Vertical boundary condition kernel (and corners)
        else if p.y==0 || p.y==M.y-1 { let t=$Y::<{1./(M.x as f32)},{1./(M.y as f32)}>(p.y); algebra::map(|i| (xy{x:0,y:t[i].0}, t[i].1)) } // Horizontal boundary condition kernel (+corners)
        else { Zero::zero() }
    } BCxy::<M> }
)}
#[macro_export] macro_rules! BC { ($X:ident, $Y:ident) => ( $crate::mesh::Rows::new($crate::BCxy!($X,$Y)) )  }
#[macro_export] macro_rules! rows_BC { ($X:ident, $Y:ident) => ( $crate::mesh::rows::<_,M,3>($crate::BCxy!($X,$Y)) ) }
