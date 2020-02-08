pub type Len = usize;
pub type Idx = usize;
use framework::core::Zero;
#[derive(Clone)] pub struct Array<T=f32, const N:Len>(pub [T; N]);
impl<T:std::fmt::Debug,const N:Len> std::fmt::Debug for Array<T,N> { fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { f.debug_list().entries(self.iter()).finish() } }
impl<T,const N:Len> std::ops::Deref for Array<T,N> { type Target = [T; N]; fn deref(&self) -> &Self::Target { &self.0 } }
impl<T,const N:Len> std::ops::DerefMut for Array<T,N> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 } }
impl<T,const N:Len> std::ops::Index<Idx> for Array<T,N> { type Output=T; fn index(&self, i:Idx) -> &Self::Output { &self.0[i] } }
impl<T,const N:Len> std::ops::IndexMut<Idx> for Array<T,N> { fn index_mut(&mut self, i:Idx) -> &mut Self::Output { &mut self.0[i] } }
impl<T:Zero, const N:Len> Zero for Array<T,N> { fn zero() -> Self { map(|_|Zero::zero()) } }
pub fn map<T, F:Fn(Idx)->T, const N:Len>(f : F) -> Array<T, N> { Array(framework::core::array::map(f)) }
impl<T:Copy,const N:Len> FnOnce<(Idx,)> for Array<T,N> { type Output=T; extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> FnMut<(Idx,)> for Array<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> Fn<(Idx,)> for Array<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self[args.0] } }

use std::ops::Mul;

pub trait Rows<const R:Len, const C:usize> = Fn(Idx)->Row<C>;
pub trait Sum : std::iter::Sum<<f32 as Mul<Self>>::Output> where f32:Mul<Self> {}
impl<T:std::iter::Sum<<f32 as Mul<Self>>::Output>> Sum for T where f32:Mul<Self> {}
pub type Row<const C:usize> = Array<(u32,f32), C>;
pub fn dot<T:Copy+Sum,const R:usize, const C:Len>(a:&Row<C>, b:&Array<T,R>) -> T where f32:Mul<T> { a.iter().map(|(j, v)| *v*b[*j as usize]).sum() }
pub trait Vector<T=f32,const N:Len> = Fn(Idx)->T;
pub fn mul<'a, T:Copy+Sum, F:Rows<R,C>, const R:usize, const C:Len>(a:&'a F, b:&'a Array<T,R>) -> impl Vector<T,R> + 'a where f32:Mul<T> {
    move |i|->T { dot(&a(i),b) }
}
/*impl<T:Copy+Sum,const R:usize, const C:Len> Mul<&Array<T,R>> for &Rows<R,C> where f32:Mul<T> {
    type Output = Array<T,R>;
    fn mul(self, b:&Array<T,R>) -> Self::Output { map(mul(self, b)) }
}*/

pub type Column<const R:usize> = Array<(u32,f32), R>;
pub trait Columns<const R:usize, const C:Len> = Fn(Idx)->Column<R>;
pub fn columns<F:Rows<K,N>, const K:usize, const N:Len, const R:usize>(rows:F) -> impl Columns<R,N> {
    move |j| {
        let mut c : Column<R> = Zero::zero(); let mut r = 0;
        for i in 0..N { for (ki, v) in rows(j).iter() { if *ki==i as u32 { c[r] = (i as u32, *v); r+=1; break; } } }
        c
    }
}

mod umfpack {
    #![allow(dead_code,non_camel_case_types,non_upper_case_globals,improper_ctypes)]
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

pub struct LU {
    column_pointers : Vec<i32>,
    row_indices : Vec<i32>,
    values : Vec<f64>,
    numeric : *mut std::os::raw::c_void // fixme: moveonly+drop
}
impl LU {
    #[allow(non_snake_case)]
    pub fn new<M:Fn(Idx)->Column<R>,const R:usize>(C:usize, A : M) -> Self {
        let mut column_pointers : Vec<i32> = vec!(0;C+1);
        let mut row_indices : Vec<i32> = Vec::with_capacity(C*R);
        let mut values : Vec<f64> = Vec::with_capacity(C*R);
        let time = std::time::Instant::now();
        for j in 0..C {
            column_pointers[j] = values.len() as i32;
            for (i, v) in A(j).iter() {
                if *v != 0. {
                    row_indices.push(*i as i32);
                    values.push(*v as f64); // fixme: might as well generate f64
                }
            }
        }
        column_pointers[C] = values.len() as i32;
        //framework::log!("eval", R, C, time.elapsed().as_millis(), values.len(), values.capacity()); let time = std::time::Instant::now();
        let mut symbolic : *mut std::os::raw::c_void = std::ptr::null_mut();
        unsafe{umfpack::umfpack_di_symbolic(C as i32, C as i32, column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), &mut symbolic, std::ptr::null(), std::ptr::null_mut())};
        //framework::log!("symbolic", time.elapsed().as_millis()); let time = std::time::Instant::now();
        let mut numeric : *mut std::os::raw::c_void = std::ptr::null_mut();
        unsafe{umfpack::umfpack_di_numeric(column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), symbolic, &mut numeric, std::ptr::null(), std::ptr::null_mut())};
        //framework::log!("numeric", time.elapsed().as_millis());
        framework::log!("LU", time.elapsed().as_millis());
        Self{column_pointers, row_indices, values, numeric}
    }
}
pub trait Solve<T> {
    fn solve<B:Vector<T,N>, const N:Len>(&self, b:B) -> Array<T,N>;
}
impl Solve<f32> for LU {
    fn solve<B:Vector<f32,N>, const N:Len>(&self, b:B) -> Array<f32,N> {
        let b : [f64;N] = framework::core::array::map(|i| b(i) as f64);
        let mut x : [f64;N] = Zero::zero(); // fixme: uninitialized
        unsafe{umfpack::umfpack_di_solve(umfpack::UMFPACK_A as i32, self.column_pointers.as_ptr(), self.row_indices.as_ptr(), self.values.as_ptr(), x.as_mut_ptr(), b.as_ptr(), self.numeric, std::ptr::null(), std::ptr::null_mut())};
        Array(framework::core::array::map(|i| x[i] as f32))
    }
}

mod tuple_vector {
#[allow(non_camel_case_types)] #[derive(Clone, Copy, Debug, PartialEq, Eq)] pub struct v2<T>(pub T, pub T);
impl<T:Copy> From<T> for v2<T> { fn from(v: T) -> Self { v2(v,v) } }
use {std::ops::{Add,Sub,Mul}, framework::core::Zero};
impl<T:Copy+Zero> Zero for v2<T> { fn zero() -> Self { T::zero().into() } }
impl<T:Add> Add<v2<T>> for v2<T> { type Output=v2<T::Output>; fn add(self, b: v2<T>) -> Self::Output { v2(self.0+b.0, self.1+b.1) } }
impl<T:Sub> Sub<v2<T>> for v2<T> { type Output=v2<T::Output>; fn sub(self, b: v2<T>) -> Self::Output { v2(self.0-b.0, self.1-b.1) } }
impl<T:std::ops::AddAssign> std::ops::AddAssign<v2<T>> for v2<T> { fn add_assign(&mut self, b: v2<T>) { self.0+=b.0; self.1+=b.1 } }
impl<T:Copy+Zero+Add<Output=T>> std::iter::Sum<v2<T>> for v2<T> { fn sum<I:Iterator<Item=v2<T>>>(iter: I) -> Self { iter.fold(Zero::zero(), Add::add) } }
fn mul<T:Copy+Mul>(a: T, b: v2<T>) -> v2<T::Output> { v2(a*b.0, a*b.1) }
impl Mul<v2<f32>> for f32 { type Output=v2<f32>; fn mul(self, b: v2<f32>) -> Self::Output { mul(self, b) } }
}
pub use tuple_vector::v2;

impl Solve<v2<f32>> for LU {
    fn solve<B:Vector<v2<f32>,N>, const N:Len>(&self, b:B) -> Array<v2<f32>,N> {
        let x0 : Array<f32, N> = self.solve(|i| b(i).0);
        let x1 : Array<f32, N> = self.solve(|i| b(i).1);
        Array(framework::core::array::map(|i| v2(x0[i], x1[i])))
    }
}
