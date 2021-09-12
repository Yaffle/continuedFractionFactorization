/*jshint esversion: 11*/
/*global BigInt*/
"use strict";

function ngcd(a, b) {
  while (b > 0) {
    const r = a - Math.floor(a / b) * b;
    a = b;
    b = r;
  }
  return a;
}

const MAX_SAFE_INTEGER = BigInt(Number.MAX_SAFE_INTEGER);
function gcd(a, b) {
  while (b > MAX_SAFE_INTEGER || a > MAX_SAFE_INTEGER && b > 0n) {
    const r = a % b;
    a = b;
    b = r;
  }
  if (b > 0n) {
    return BigInt(ngcd(Number(a), Number(b)));
  }
  return a;
}

function log2(x) {
  return BigInt(x.toString(16).length * 4);
}
function sqrt(x) {
  if (x < 2n**52n) {
    return BigInt(Math.floor(Math.sqrt(Number(x) + 0.5)));
  }
  const q = (log2(x) / 4n);
  const initialGuess = ((sqrt(x >> (q * 2n)) + 1n) << q);
  var a = initialGuess, b = a+1n;
  while(a<b) {
    b = a;
    a = (b + x/b)/2n;
  }
  return b;
}

// See https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html

function* continuedFractionForSqrt(n) {
  function floorDiv(a, b) {
    return typeof a === "bigint" && typeof b === "bigint" ? a / b : Math.floor(a / b);
  }
  n = BigInt(n);
  if (n < 0n) {
    throw new RangeError();
  }
  const root = BigInt(sqrt(n));
  const remainder = n - root * root;
  const i = Number(root) * 2 <= Number.MAX_SAFE_INTEGER ? Number(root) : root;
  yield i;
  if (remainder !== 0n) {
    // sqrt(n) = floor(sqrt(n)) + 1 / x
    // expand continued fraction:
    // x_k = floor(x_k) + 1 / x_(k+1)
    // each x_k will have form x_k = (sqrt(n) + y_k) / z_k
    const one = i / i;
    let zprev = one;
    let z = typeof i === "number" ? Number(remainder) : remainder;
    let y = i;
    do {
      const q = floorDiv((i + y), z);
      const ynew = q * z - y;
      const znew = zprev + q * (y - ynew);
      y = ynew;
      zprev = z;
      z = znew;
      //console.assert(one <= q && q <= i + i);
      //console.assert(one <= y && y <= i);
      //console.assert(one <= z && z <= i + i);
      yield q;
    } while (zprev !== one); //TODO: why is it cycling here - ?
    console.assert(y === i && BigInt(z) === remainder && zprev === one);
  }
}

function* sqrtConvergentNumeratorsModN(N) {
  // https://en.wikipedia.org/wiki/Continued_fraction#:~:text=The%20successive%20convergents%20are%20given%20by%20the%20formula
  N = BigInt(N);
  let hprev = 0n;
  let h = 1n;
  for (const a of continuedFractionForSqrt(N)) { // TODO: why do we stop after the first cycle ?
    [hprev, h] = [h, BigInt(a) * h + hprev];
    // optimization from the https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html :
    h = h % N;
    yield h;
  }
}

function modPow(a, n, m) {
  a = Number(a);
  n = Number(n);
  m = Number(m);
  if (Math.abs(m) > Math.floor(Math.sqrt(Number.MAX_SAFE_INTEGER)) || n < 0 || a % m === 0) {
    throw new RangeError();
  }
  a = a % m;
  let result = 1;
  while (n > 0) {
    if ((n % 2) === 1) {
      result = (result * a) % m;
    }
    n = Math.floor(n / 2);
    a = (a**2) % m;
  }
  return result;
}

function isQuadraticResidueModuloPrime(a, p) {
  if (p === 2) {
    // "Modulo 2, every integer is a quadratic residue." - https://en.wikipedia.org/wiki/Quadratic_residue#Prime_modulus
    return true;
  }
  console.assert(p % 2 === 1);
  // https://en.wikipedia.org/wiki/Euler%27s_criterion
  const amodp = Number(BigInt(a) % BigInt(p));
  if (amodp === 0) {
    return true;
  }
  const value = modPow(amodp, (p - 1) / 2, p);
  console.assert(value === 1 || value === p - 1);
  return value === 1;
}

function L(N) {  // exp(sqrt(log(n)*log(log(n))))
  const e = Math.max(N.toString(16).length * 4 - 4 * 12, 0);
  const lnn = Math.log(Number(N >> BigInt(e))) + Math.log(2) * e;
  return Math.exp(Math.sqrt(lnn * Math.log(lnn)));
}

function product(array) {
  //return array.reduce((p, a) => p * BigInt(a), 1n);
  if (array.length === 0) {
    return 1n;
  }
  if (array.length === 1) {
    return BigInt(array[0]);
  }
  const m = Math.floor(array.length / 2);
  return product(array.slice(0, m)) * product(array.slice(m));
}

const T1 = (1n << 64n);
function isSmoothOverProduct(a, product, product1) {
  if (a >= T1) {
    var g1 = BigInt(gcd(a, product1));
    if (g1 > 1n) {
      a /= g1;
    }
  }
  if (a === 0n) {
    throw new RangeError();
  }
  // quite test from https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=is_smooth_over_prod
  let g = BigInt(gcd(a, product % a));
  while (g > 1n) {
    a /= g;
    g = BigInt(gcd(a, g));
  }
  return a;
}

function getSmoothFactorization(a, base) {  
  if (a === 0n) {
    return [];//TODO: ?
  }
  var value = BigInt(a);
  var result = new Array(base.length);
  for (var i = 0; i < base.length; i += 1) {
    var p = base[i];
    var e = 0;
    while (value % BigInt(p) === 0n) {
      value /= BigInt(p);
      e += 1;
    }
    result[i] = e;
  }
  return value === 1n ? result : null;
}

// (X**2 - Y) % N === 0, where Y is a smooth number
function CongruenceOfsquareOfXminusYmoduloN(X, Y, N) {
  this.X = X;
  this.Y = Y;
  this.N = N;
}
CongruenceOfsquareOfXminusYmoduloN.prototype.toString = function () {
  const X = this.X, Y = this.Y, N = this.N;
  return 'X**2 ≡ Y (mod N)'.replaceAll('X', X).replaceAll('Y', Y).replaceAll('N', N);
};

function getCongruencesUsingContinuedFraction(primes, n) {
  const primesProduct = product(primes);
  const product1 = BigInt(primes.reduce((p, prime) => p * prime <= Number.MAX_SAFE_INTEGER ? p * prime : p, 1));
  n = BigInt(n);
  const g = sqrtConvergentNumeratorsModN(n);
  let congruences = [];
  while (congruences.length < primes.length + 1) {
    const A_k = g.next().value;
    if (A_k == null) {
      return congruences;//?
    }
    const X = A_k % n; // A_k mod n
    let Y = (X * X) % n; // (A_k)^2 mod n
    if (Y > n - Y) {//TODO: ???
      Y = n - Y;
    }
    //console.log(X, Y);
    if (Y === 0n) {
      return [new CongruenceOfsquareOfXminusYmoduloN(X, Y, n)];//TODO: ???
    }
    const s = isSmoothOverProduct(Y, primesProduct, product1);
    if (s === 1n) {
      congruences.push(new CongruenceOfsquareOfXminusYmoduloN(X, Y, n));
    }
  }
  return congruences;
}

function BitSet(size) {
  this.data = new Array(Math.ceil(size / 32)).fill(0);
  this.size = size;
}
BitSet.prototype.has = function (index) {
  if (index >= this.size) {
    throw new RangeError();
  }
  return (this.data[Math.floor(index / 32)] & (1 << (index % 32))) !== 0;
};
BitSet.prototype.add = function (index) {
  if (index >= this.size) {
    throw new RangeError();
  }
  this.data[Math.floor(index / 32)] |= (1 << (index % 32));
};
BitSet.prototype.xor = function (other) {
  const n = other.data.length;
  if (n !== this.data.length) {
    throw new RangeError();
  }
  for (let i = 0; i < n; i += 1) {
    this.data[i] ^= other.data[i];
  }
};
BitSet.prototype.toString = function () {
  return this.data.map(x => (x >>> 0).toString(2).padStart(32, '0').split('').reverse().join('')).join('');
};

// factorizations - array of arrays of powers
// returns combinations of factorizations which result in zero vector by modulo 2
function solve(factorizations) {
  const primeBaseSize = factorizations.length === 0 ? 0 : factorizations[0].length;
  // 1. build the augmented matrix
  // 2. row reduce
  var M = new Array(factorizations.length);
  for (var i = 0; i < M.length; i++) {
    const row = new BitSet(primeBaseSize + M.length);
    const f = factorizations[i];
    for (var j = 0; j < primeBaseSize; j += 1) {
      if (f[j] % 2 !== 0) {
        row.add(j);
      }
    }
    row.add(primeBaseSize + i);
    M[i] = row;
  }
  //console.log(M.map(x => x.toString()).join('\n'))
  var pivotRow = 0;
  for (var pivotColumn = 0; pivotColumn < primeBaseSize; pivotColumn += 1) {
    var row = pivotRow;
    while (row < M.length && !M[row].has(pivotColumn)) {
      row += 1;
    }
    if (row < M.length) {
      if (row !== pivotRow) {
        // swap rows:
        const tmp1 = M[row];
        M[row] = M[pivotRow];
        M[pivotRow] = tmp1;
      }
      // row-reduction:
      for (var i = pivotRow + 1; i < M.length; i++) {
        if (M[i].has(pivotColumn)) {
          M[i].xor(M[pivotRow]);
        }
      }
      pivotRow += 1;
    }
  }
  var solutions = [];
  while (pivotRow < M.length) {
    var row = M[pivotRow];
    var solution = [];
    for (var i = 0; i < M.length; i += 1) {
      if (row.has(primeBaseSize + i)) {
        solution.push(i);
      }
    }
    solutions.push(solution);
    pivotRow += 1;
  }
  return solutions;
}

function primes(MAX) {
  let sieve = new Array(MAX + 1).fill(true);
  let result = [];
  result.push(2);
  for (let i = 3; i <= MAX; i += 2) {
    if (sieve[i]) {
      result.push(i);
      for (let j = i * i; j <= MAX; j += 2 * i) {
        sieve[j] = false;
      }
    }
  }
  return result;
}

function ContinuedFractionFactorization(N) {
  N = BigInt(N);
  //if (isPrime(N)) {
  //  return N;
  //}
  for (let k = 1n;; k += 1n) {
    const kN = k * N;
    // https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=optimal%20value :
    const B = Math.min(Math.floor(Math.sqrt(L(kN))), 2**25 - 1);
    // as we are check for smoothness the number P_k^2 - N * Q_k^2, where P_k/Q_k is a convergent,
    // if p is a factor of this number then P_k^2 - N * Q_k^2 = 0 (mod p),
    // after multiplication by Q_k^-2 we have (Q_k^-1*P_k)^2 = N (mod p),
    // x^2 = N (mod p), so N is a quadratic residue modulo N
    // см. "Теоретико-числовые методы в криптографии" (О.Н. Герман, Ю.В. Нестеренко), страница 202
    const primeBase = primes(B).filter(p => isQuadraticResidueModuloPrime(kN, p));
    const congruences = getCongruencesUsingContinuedFraction(primeBase, kN); // congruences X_k^2 = Y_k mod N, where Y_k is smooth over the prime base
    const solutions = solve(congruences.map(c => getSmoothFactorization(c.Y, primeBase))); // find products of Y_k = Y, so that Y^2 is a perfect square
    for (const solution of solutions) {
      const X = product(solution.map(i => BigInt(congruences[i].X)));
      const Y = product(solution.map(i => BigInt(congruences[i].Y))); // = sqrt(X**2 % N)
      const x = X;
      const y = BigInt(sqrt(Y));
      console.assert(y**2n === BigInt(Y));
      const g = BigInt(gcd(x + y, N));
      if (g !== 1n && g !== N) {
        return g;
      }
    }
  }
}

ContinuedFractionFactorization.testables = {
  continuedFractionForSqrt,
  sqrtConvergentNumeratorsModN,
  primes: primes,
  L: L,
  getSmoothFactorization: getSmoothFactorization,
  CongruenceOfsquareOfXminusYmoduloN: CongruenceOfsquareOfXminusYmoduloN,
  solve: solve,
  product: product,
  isQuadraticResidueModuloPrime: isQuadraticResidueModuloPrime
};

export default ContinuedFractionFactorization;
