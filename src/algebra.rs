pub type Len = usize;
pub type Idx = usize;
use framework::core::Zero;
#[derive(Clone)] pub struct Array<T=f32, const N:Len>(pub [T; N]);
impl<T,const N:Len> std::ops::Deref for Array<T,N> { type Target = [T; N]; fn deref(&self) -> &Self::Target { &self.0 } }
impl<T,const N:Len> std::ops::DerefMut for Array<T,N> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 } }
impl<T,const N:Len> std::ops::Index<Idx> for Array<T,N> { type Output=T; fn index(&self, i:Idx) -> &Self::Output { &self.0[i] } }
impl<T,const N:Len> std::ops::IndexMut<Idx> for Array<T,N> { fn index_mut(&mut self, i:Idx) -> &mut Self::Output { &mut self.0[i] } }
pub fn map<T, F:Fn(Idx)->T, const N:Len>(f : F) -> Array<T, N> { Array(framework::core::array::map(f)) }
impl<T:Copy,const N:Len> FnOnce<(Idx,)> for Array<T,N> { type Output=T; extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> FnMut<(Idx,)> for Array<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> Fn<(Idx,)> for Array<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self[args.0] } }
impl<T:Zero, const N:Len> Zero for Array<T,N> { fn zero() -> Self { map(|_|Zero::zero()) } }

use std::ops::Mul;

// Vector
pub trait Vector<T=f32,const N:Len> = Fn(Idx)->T;
pub type BoxVector<'a,T,const N:Len> = crate::compose::BoxFn<'a,(Idx,),T>;
// impl scalar · vector for Array
//pub fn mul<'a, T:Copy, const N:Len>(a:f32, b:&'a Array<T,N>) -> BoxVector<'a,T,N> where f32:Mul<T,Output=T> { BoxVector::new(move |i| a*b[i]) }

/// Operator
pub trait Operator<T=f32, const N:Len> = Fn(&dyn Vector<T,N>)->BoxVector<T,N>;
pub trait OperatorOnce<T=f32, const N:Len> = FnOnce(&dyn Vector<T,N>)->BoxVector<T,N>;
pub type BoxOperator<'a, T=f32, const N:Len> = crate::compose::BoxFn<'a, (&'a dyn Vector<T,N>,), BoxVector<'a,T,N>>;

// vector · operator
impl<'a, T, const N:Len> Mul<&'a dyn Operator<T,N>> for &'a Array<f32,N> where f32:Mul<T,Output=T> {
    type Output = BoxOperator<'a,T,N>;
    fn mul(self, b:&'a dyn Operator<T,N>) -> Self::Output { BoxOperator::new(move |v| { let b=b(v); BoxVector::new(move |i| self[i] * b(i)) }) }
}

pub trait Dot<T> { fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T; }
//pub trait Sum : std::iter::Sum<<f32 as Mul<Self>>::Output> where f32:Mul<Self> {}
pub trait Sum : std::iter::Sum<<f32 as Mul<Self>>::Output>+Zero+std::ops::AddAssign<<f32 as std::ops::Mul<Self>>::Output> where f32:Mul<Self> {}
impl<T:std::iter::Sum<<f32 as Mul<T>>::Output>+Zero+std::ops::AddAssign<<f32 as std::ops::Mul<T>>::Output>> Sum for T where f32:Mul<T> {}
impl<T:Sum> Dot<T> for std::ops::Range<Idx> where f32:Mul<T> {
    fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T {
        //self.map(|k| a(k) * b(k)).sum() // ICE `Unimplemented` selecting `Binder...`
        let mut sum = Zero::zero();
        for k in self { sum += a(k) * b(k) }
        sum
    }
}
impl<T:Sum> Dot<T> for std::ops::RangeTo<Idx> where f32:Mul<T> {
    fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T { (0..self.end).dot(a,b) }
}
/*pub trait DotN<T> { fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized,const N:Len>(self, a : &A, b : &B) -> T; }
impl<T> DotN<T> for std::ops::RangeFrom<Idx> where std::ops::Range<Idx>:Dot<T> {
    fn dot<A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(self, a : &A, b : &B) -> T { (self.start..N).dot(a,b) }
}
impl<T> DotN<T> for std::ops::RangeFull where std::ops::RangeFrom<Idx>:DotN<T> {
    fn dot<A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(self, a : &A, b : &B) -> T { (0..).dot(a,b) }
}*/
//pub fn dot<T,A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(a : &A, b : &B) -> T where std::ops::RangeFull:DotN<T> { (..).dot(a,b) }
pub fn dot<T:Sum,A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(a : &A, b : &B) -> T where std::ops::RangeTo<Idx>:Dot<T>, f32:Mul<T> { (..N).dot(a,b) }

pub trait Matrix<const N:Len> = Fn(Idx,Idx)->f32;
pub fn matrix_mul<'a, A:Matrix<N>+'a,T:Sum,B:Vector<T,N>+?Sized,const N:Len>(a:A, b:&'a B) -> impl Vector<T,N> + 'a where f32:Mul<T> {
    move |i| dot(&|j| a(i, j), b)
}

/*pub struct SparseVector<T=f32,const N:Len> {
    indices : Vec<Idx>,
    values : Vec<T>
}
// derive[Default] req T:Default for no reason
impl<T,const N:Len> Default for SparseVector<T,N> { fn default() -> Self { Self{indices:Default::default(),values:Default::default()} } }

impl<T,const N:Len> SparseVector<T,N> {
    fn get(&self, i : Idx) -> Option<T> where T:Copy { Some( self.values[ self.indices.binary_search(&i).ok()? ] ) }
    fn insert(&mut self, i : Idx, v : T) -> Option<T> {
        assert!(i<N);
        match self.indices.binary_search(&i) {
            Ok(index) => Some(std::mem::replace(&mut self.values[index], v)),
            Err(index) => { self.indices.insert(index, i); self.values.insert(index, v); None }
        }
    }
}

impl<T:Copy+Zero,const N:Len> FnOnce<(Idx,)> for SparseVector<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> FnMut<(Idx,)> for SparseVector<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> Fn<(Idx,)> for SparseVector<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self.get(args.0).unwrap_or(T::zero()) } }

pub struct CSC<T=f32,const N:Len> /*where [SparseVector<T,N>;N]:*/ { columns : [SparseVector<T,N>; N] }
impl<T,const N:Len> Zero for CSC<T,N> { fn zero() -> Self { Self{columns: framework::core::array::map(|_|Default::default())} } }

impl<T:Copy+Zero,const N:Len> CSC<T,N> {
    fn column(&self, j : Idx) -> &SparseVector<T,N> { &self.columns[j as usize] }
    fn column_mut(&mut self, j : Idx) -> &mut SparseVector<T,N> { &mut self.columns[j as usize] }
    fn get(&self, i : Idx, j : Idx) -> Option<T> { self.column(j).get(i) }
    fn value(&self, i : Idx, j : Idx) -> T { self.get(i, j).unwrap_or(T::zero()) }
    fn line(&self, i : Idx) -> impl Vector<T,N> + '_ { move |j| self.value(i, j) }
    fn insert(&mut self, i:Idx, j:Idx, v:T) -> Option<T> { self.column_mut(j).insert(i, v) }
    fn insert_once(&mut self, i:Idx, j:Idx, v:T) where T:std::fmt::Debug { self.insert(i, j, v).unwrap_none() }
}

impl<T:Copy+Zero,const N:Len> FnOnce<(Idx,Idx)> for CSC<T,N> { type Output=T; extern "rust-call" fn call_once(self, args: (Idx,Idx)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> FnMut<(Idx,Idx)> for CSC<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,Idx)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> Fn<(Idx,Idx)> for CSC<T,N> { extern "rust-call" fn call(&self, args: (Idx,Idx)) -> Self::Output { self.value(args.0, args.1)} }

pub type LU<const N:Len> = CSC<f32,N>;
impl<const N:Len> LU<N> {
  #[allow(non_snake_case)]
    pub fn new<M:Matrix<N>>(A : M) -> Box<Self> {
        let time = std::time::Instant::now();
        let mut LU : Box<Self> = unsafe { Box::new_zeroed().assume_init() };
        for i in 0..N {
            for j in i..N { let v = A(i, j) - (..i).dot(&LU.line(i), LU.column(j)); LU.insert_once(i, j, v);  }
            for j in i+1..N { let v = (1. / LU.get(i, i).unwrap()) * (A(j, i) - (..i).dot(&LU.line(i), LU.column(j))); LU.insert_once(j, i, v); }
        }
        LU
    }
    pub fn solve<T:Copy+Sum+std::ops::Sub<Output=T>,B:Vector<T,N_>, const N_:Len>(&self, b:B) -> Array<T,N_> where f32:Mul<T,Output=T> {
        assert_eq!(N, N_);
        let mut y : Array<T,N> = Zero::zero(); // Ly = b (uninit ok)
        for i in 0..N { y[i] = b(i) - (..i).dot(&self.line(i), &y); }
        let mut x : Array<T,N_> = Zero::zero(); // Ux = y (uninit ok)
        for i in N-1..=0 { x[i] = (1. / self(i, i)) * (y[i] - (i+1..N).dot(&self.line(i), &x)); }
        x
    }
}*/

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
    //pub fn new<M:Matrix<N>, const N:Len>(A : M) -> Self {
    pub fn new<M:Fn(Idx,Idx)->f32>(N:Len, A : M) -> Self {
        let time = std::time::Instant::now();
        //let mut column_pointers : [i32; N+1] = Zero::zero();
        let mut column_pointers : Vec<i32> = Vec::new(); for _ in 0..N+1 { column_pointers.push(0); }
        let mut row_indices : Vec<i32> = Vec::new();
        let mut values : Vec<f64> = Vec::new();
        for j in 0..N {
            column_pointers[j] = values.len() as i32;
            for i in 0..N {
                let v = A(i,j);
                if v != 0. {
                    row_indices.push(i as i32);
                    values.push(v as f64); // fixme: might as well generate f64
                }
            }
        }
        column_pointers[N] = values.len() as i32;
        let mut symbolic : *mut std::os::raw::c_void = std::ptr::null_mut();
        unsafe{umfpack::umfpack_di_symbolic(N as i32, N as i32, column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), &mut symbolic, std::ptr::null(), std::ptr::null_mut())};
        let mut numeric : *mut std::os::raw::c_void = std::ptr::null_mut();
        unsafe{umfpack::umfpack_di_numeric(column_pointers.as_ptr(), row_indices.as_ptr(), values.as_ptr(), symbolic, &mut numeric, std::ptr::null(), std::ptr::null_mut());
        framework::log!(time.elapsed().as_millis())};
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

#[allow(non_camel_case_types)] #[derive(Clone, Copy, Debug, PartialEq, Eq)] pub struct v2<T>(pub T, pub T);
impl<T:Copy> From<T> for v2<T> { fn from(v: T) -> Self { v2(v,v) } }
impl<T:Copy+Zero> Zero for v2<T> { fn zero() -> Self { T::zero().into() } }
use std::ops::{Add,Sub};
impl<T:Add> Add<v2<T>> for v2<T> { type Output=v2<T::Output>; fn add(self, b: v2<T>) -> Self::Output { v2(self.0+b.0, self.1+b.1) } }
impl<T:Sub> Sub<v2<T>> for v2<T> { type Output=v2<T::Output>; fn sub(self, b: v2<T>) -> Self::Output { v2(self.0-b.0, self.1-b.1) } }
impl<T:std::ops::AddAssign> std::ops::AddAssign<v2<T>> for v2<T> { fn add_assign(&mut self, b: v2<T>) { self.0+=b.0; self.1+=b.1 } }
impl<T:Copy+Zero+Add<Output=T>> std::iter::Sum<v2<T>> for v2<T> { fn sum<I:Iterator<Item=v2<T>>>(iter: I) -> Self { iter.fold(Zero::zero(), Add::add) } }
fn mul<T:Copy+Mul>(a: T, b: v2<T>) -> v2<T::Output> { v2(a*b.0, a*b.1) }
impl Mul<v2<f32>> for f32 { type Output=v2<f32>; fn mul(self, b: v2<f32>) -> Self::Output { mul(self, b) } }

//#[allow(non_camel_case_types)] pub type v2xf32 = v2<f32>;

impl Solve<v2<f32>> for LU {
    fn solve<B:Vector<v2<f32>,N>, const N:Len>(&self, b:B) -> Array<v2<f32>,N> {
        let x0 : Array<f32, N> = self.solve(|i| b(i).0);
        let x1 : Array<f32, N> = self.solve(|i| b(i).1);
        Array(framework::core::array::map(|i| v2(x0[i], x1[i])))
    }
}
