use framework::{assert,core::Zero};
#[derive(Clone)] pub struct Array<T=f32, const N:usize>(pub [T; N]);
impl<T:std::fmt::Debug,const N:usize> std::fmt::Debug for Array<T,N> { fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { f.debug_list().entries(self.iter()).finish() } }
impl<T,const N:usize> std::ops::Deref for Array<T,N> { type Target = [T; N]; fn deref(&self) -> &Self::Target { &self.0 } }
impl<T,const N:usize> std::ops::DerefMut for Array<T,N> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 } }
/*impl<T,I,const N:usize> std::ops::Index<I> for Array<T,N> where [T;N]:std::ops::Index<I> {
    type Output=<[T;N] as std::ops::Index<I>>::Output; fn index(&self, i:I) -> &Self::Output { &self.0[i] }
}
impl<T,I,const N:usize> std::ops::IndexMut<I> for Array<T,N> where [T;N]:std::ops::IndexMut<I> { fn index_mut(&mut self, i:I) -> &mut Self::Output { &mut self.0[i] } }*/
impl<T,const N:usize> std::ops::Index<usize> for Array<T,N> { type Output=T; fn index(&self, i:usize) -> &Self::Output { &self.0[i] } }
impl<T,const N:usize> std::ops::IndexMut<usize> for Array<T,N> { fn index_mut(&mut self, i:usize) -> &mut Self::Output { &mut self.0[i] } }
impl<T,const N:usize> std::ops::Index<std::ops::Range<usize>> for Array<T,N> { type Output=[T]; fn index(&self, i:std::ops::Range<usize>) -> &Self::Output { &self.0[i] } }
impl<T,const N:usize> std::ops::IndexMut<std::ops::Range<usize>> for Array<T,N> { fn index_mut(&mut self, i:std::ops::Range<usize>) -> &mut Self::Output { &mut self.0[i] } }

impl<T:Zero, const N:usize> Zero for Array<T,N> { fn zero() -> Self { map(|_|Zero::zero()) } }
pub type Idx = usize;
pub fn map<T, F:Fn(Idx)->T, const N:usize>(f : F) -> Array<T, N> { Array(framework::core::array::map(f)) }
impl<T:Copy,const N:usize> FnOnce<(Idx,)> for Array<T,N> { type Output=T; extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:usize> FnMut<(Idx,)> for Array<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:usize> Fn<(Idx,)> for Array<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self[args.0] } }

use std::ops::Mul;

pub trait Rows<const R:usize, const C:usize> = Fn(Idx)->Row<C>;
pub trait Sum : std::iter::Sum<<f32 as Mul<Self>>::Output> where f32:Mul<Self> {}
impl<T:std::iter::Sum<<f32 as Mul<Self>>::Output>> Sum for T where f32:Mul<Self> {}
pub type Row<const C:usize> = Array<(u32,f32), C>;
pub fn dot<T:Copy+Sum,const R:usize, const C:usize>(a:&Row<C>, b:&Array<T,R>) -> T where f32:Mul<T> { a.iter().map(|(j, v)| *v*b[*j as usize]).sum() }
pub trait Vector<T=f32,const N:usize> = Fn(Idx)->T;
pub fn mul<'t, T:Copy+Sum, F:Rows<R,C>, const R:usize, const C:usize>(a:&'t F, b:&'t Array<T,R>) -> impl Vector<T,R>+'t where f32:Mul<T> {
    move |i|->T { dot(&a(i),b) }
}
/*impl<T:Copy+Sum,const R:usize, const C:usize> Mul<&Array<T,R>> for &Rows<R,C> where f32:Mul<T> {
    type Output = Array<T,R>;
    fn mul(self, b:&Array<T,R>) -> Self::Output { map(mul(self, b)) }
}*/

pub type Column<const R:usize> = Array<(u32,f32), R>;
//pub trait Columns<const R:usize, const C:usize> = Fn(Idx)->Column<R>;
pub type Columns<const R:usize, const C:usize> = Array<Column<R>, C>;
pub fn columns<F:Rows<N,K>, const K:usize, const N:usize, const R:usize>(rows:F) -> /*impl*/ Box<Columns<R,N>> {
    /*move |j:usize| {
        let mut c : Column<R> = Zero::zero(); let mut r = 0;
        for i in 0..N {
            for (kj, v) in rows(i).iter() {
                if *kj==j as u32 {
                    /*let mut k = 0; while c[k].0 < i { k+=1; }
                    //c[k+1..r+1] = c[k..r];
                    impl<T,const N:usize> Array<T,N> {
                        unsafe fn insert(&mut self, len: usize, index: usize, element: T) {
                            let p = self.as_mut_ptr().add(index);
                            std::ptr::copy(p, p.offset(1), len - index);
                            std::ptr::write(p, element);
                        }
                    }
                    //c[r] = (i, *v);
                    unsafe{c.insert(r, k, (i, *v))};*/
                    c[r] = (i as u32, *v);
                    r+=1;
                    break;
                }
            }
        }
        //for i in 0..c.len() { if c[i].1 != 0. { for j in 0..i { assert!(c[j].0 < c[i].0, c); } } }
        c
    }*/
    let mut columns : Box<Array<Column<R>, N>> = box Zero::zero(); // fixme: may be uninitialized
    let mut column_lens = vec!(0;N); //[0;N];
    for i in 0..N {
        for (j, v) in rows(i).iter() {
            if *v != 0. {
                let len = &mut column_lens[*j as usize];
                assert!(*len < R, columns[*j as usize]);
                columns[*j as usize][*len as usize] = (i as u32, *v);
                *len+=1;
            }
        }
    }
    //move |j:usize| { columns[j].clone() }
    columns
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
    //pub fn new<M:Fn(Idx)->Column<R>,const R:usize>(C:usize, A : M) -> Self {
    pub fn new<const C:usize, const R:usize>(A : &Columns<R,C>) -> Self {
        let mut column_pointers : Vec<i32> = vec!(0;C+1);
        let mut row_indices : Vec<i32> = Vec::with_capacity(C*R);
        let mut values : Vec<f64> = Vec::with_capacity(C*R);
        for j in 0..C {
            column_pointers[j] = values.len() as i32;
            //for (i, v) in A(j).iter() {
            for (i, v) in A[j].iter() {
                //std::assert!(*v!=0.);
                if *v != 0. {
                    //assert!(!row_indices[column_pointers[j] as usize..].contains(&(*i as i32)), &row_indices[column_pointers[j] as usize..], A(j));
                    //for ri in &row_indices[column_pointers[j] as usize..] { assert!(ri < &(*i as i32), &row_indices[column_pointers[j] as usize..], A(j)); }
                    for ri in &row_indices[column_pointers[j] as usize..] { assert!(ri < &(*i as i32), &row_indices[column_pointers[j] as usize..], A[j]); }
                    row_indices.push(*i as i32);
                    values.push(*v as f64); // fixme: might as well generate f64
                }
            }
        }
        column_pointers[C] = values.len() as i32;
        let mut symbolic : *mut std::os::raw::c_void = std::ptr::null_mut();
        let status = unsafe{umfpack::umfpack_di_symbolic(C as i32, C as i32, column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), &mut symbolic, std::ptr::null(), std::ptr::null_mut())};
        assert!(status == 0, "symbolic", status);
        let mut numeric : *mut std::os::raw::c_void = std::ptr::null_mut();
        let status = unsafe{umfpack::umfpack_di_numeric(column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), symbolic, &mut numeric, std::ptr::null(), std::ptr::null_mut())};
        assert!(status == 0, "numeric", status);
        Self{column_pointers, row_indices, values, numeric}
    }
}
pub trait Solve<T> {
    fn solve<B:Vector<T,N>, const N:usize>(&self, b:B) -> Array<T,N>;
}
impl Solve<f32> for LU {
    fn solve<B:Vector<f32,N>, const N:usize>(&self, b:B) -> Array<f32,N> {
        let b : [f64;N] = framework::core::array::map(|i| b(i) as f64);
        let mut x : [f64;N] = Zero::zero(); // fixme: uninitialized
        let status = unsafe{umfpack::umfpack_di_solve(umfpack::UMFPACK_A as i32, self.column_pointers.as_ptr(), self.row_indices.as_ptr(), self.values.as_ptr(), x.as_mut_ptr(), b.as_ptr(), self.numeric, std::ptr::null(), std::ptr::null_mut())};
        assert!(status == 0, "solve", status);
        map(|i| x[i] as f32)
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
    fn solve<B:Vector<v2<f32>,N>, const N:usize>(&self, b:B) -> Array<v2<f32>,N> {
        let x0 : Array<f32, N> = self.solve(|i| b(i).0);
        let x1 : Array<f32, N> = self.solve(|i| b(i).1);
        Array(framework::core::array::map(|i| v2(x0[i], x1[i])))
    }
}
