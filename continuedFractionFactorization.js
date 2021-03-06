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

function gcd(a, b) {
  while (BigInt.asUintN(53, b) !== b) {
    const r = a % b;
    a = b;
    b = r;
  }
  if (b !== 0n) {
    if (BigInt.asUintN(53, a) !== a) {
      const r = a % b;
      a = b;
      b = r;
    }
    if (b !== 0n) {
      return BigInt(ngcd(Number(a), Number(b)));
    }
  }
  return a;
}

function log2(x) {
  return BigInt(x.toString(16).length * 4);
}

function sqrt(x) {
  if (x < BigInt((Number.MAX_SAFE_INTEGER + 1) / 2)) {
    return BigInt(Math.floor(Math.sqrt(Number(x) + 0.5)));
  }
  const q = (log2(x) / 4n);
  const initialGuess = ((sqrt(x >> (q * 2n)) + 1n) << q);
  let a = initialGuess, b = a+1n;
  while(a<b) {
    b = a;
    a = (b + x/b)/2n;
  }
  return b;
}

// See https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html

function continuedFractionForSqrt(n) {
  function floorDiv(a, b) {
    return typeof a === "number" && typeof b === "number" ? Math.floor(a / b) : a / b;
  }
  n = BigInt(n);
  if (n < 0n) {
    throw new RangeError();
  }
  const root = BigInt(sqrt(n));
  const remainder = n - root * root;
  const i = Number(root) * 2 <= Number.MAX_SAFE_INTEGER ? Number(root) : root;
  let state = 0;
  // sqrt(n) = floor(sqrt(n)) + 1 / x
  // expand continued fraction:
  // x_k = floor(x_k) + 1 / x_(k+1)
  // each x_k will have form x_k = (sqrt(n) + y_k) / z_k
  const one = i / i;
  let zprev = one;
  let z = typeof i === "number" ? Number(remainder) : remainder;
  let y = i;
  const iterator = {
    next: function continuedFractionForSqrt() {
      if (state === 0) {
        state = 1;
        return {value: i, done: false};
      }
      if (remainder !== 0n) {
        while (state === 1 || zprev !== one) { //TODO: why is it cycling here - ?
          state = 2;
          const q = floorDiv((i + y), z);
          const ynew = q * z - y;
          const znew = zprev + q * (y - ynew);
          y = ynew;
          zprev = z;
          z = znew;
          //console.assert(one <= q && q <= i + i);
          //console.assert(one <= y && y <= i);
          //console.assert(one <= z && z <= i + i);
          return {value: q, done: false};
        }
        //console.assert(y === i && BigInt(z) === remainder && zprev === one);
      }
      return {value: undefined, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

function modPowSmall(base, exponent, modulus) {
  base = Number(base);
  exponent = Number(exponent);
  modulus = Number(modulus);
  if (Math.max(Math.pow(modulus, 2), Math.pow(base, 2)) > Number.MAX_SAFE_INTEGER) {
    throw new RangeError();
  }
  let accumulator = 1;
  while (exponent !== 0) {
    if (exponent % 2 === 0) {
      exponent /= 2;
      base = (base * base) % modulus;
    } else {
      exponent -= 1;
      accumulator = (accumulator * base) % modulus;
    }
  }
  return accumulator;
}

function isQuadraticResidueModuloPrime(a, p) {
  a = BigInt(a);
  p = Number(p);
  if (p === 2) {
    // "Modulo 2, every integer is a quadratic residue." - https://en.wikipedia.org/wiki/Quadratic_residue#Prime_modulus
    return true;
  }
  // https://en.wikipedia.org/wiki/Euler%27s_criterion
  const amodp = Number(BigInt(a) % BigInt(p));
  if (amodp === 0) {
    return true;
  }
  console.assert(p % 2 === 1);
  const value = modPowSmall(amodp, (p - 1) / 2, p);
  console.assert(value === 1 || value === p - 1);
  return value === 1;
}

function log(N) {
  const e = Math.max(N.toString(16).length * 4 - 4 * 12, 0);
  const lnn = Math.log(Number(N >> BigInt(e))) + Math.log(2) * e;
  return lnn;
}

function L(N) {  // exp(sqrt(log(n)*log(log(n))))
  const lnn = log(N);
  return Math.exp(Math.sqrt(lnn * Math.log(lnn)));
}

function product(array) {
  const n = array.length;
  const m = Math.floor(n / 2);
  return n === 0 ? 1n : (n === 1 ? BigInt(array[0]) : product(array.slice(0, m)) * product(array.slice(m)));
}

function isSmoothOverProduct(a, product, product1) {
  a = a < 0n ? -a : a;
  if (BigInt.asUintN(64, a) !== a) {
    const g1 = gcd(a, product1);
    if (g1 !== 1n) {
      a /= g1;
    }
  }
  if (a === 0n) {
    throw new RangeError();
  }
  // quite test from https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=is_smooth_over_prod
  let g = gcd(a, product % a);
  while (g !== 1n) {
    a /= g;
    g = gcd(a, g);
  }
  return a;
}

function getSmoothFactorization(a, base) {  
  let value = BigInt(a);
  if (value === 0n) {
    return [];
  }
  let result = [];
  if (value < 0n) {
    result.push(-1);
    value = -value;
  }
  for (let i = 0; i < base.length; i += 1) {
    const p = BigInt(base[i]);
    while (value % p === 0n) {
      value /= p;
      result.push(p);
    }
  }
  return value === 1n ? result : null;
}

// (X**2 - Y) % N === 0, where Y is a smooth number
function CongruenceOfsquareOfXminusYmoduloN(X, Y, N, factorization) {
  this.X = X;
  this.Y = Y;
  this.N = N;
  this.factorization = factorization;
}
CongruenceOfsquareOfXminusYmoduloN.prototype.toString = function () {
  const X = this.X, Y = this.Y, N = this.N;
  return 'X**2 ??? Y (mod N), Y = F'.replaceAll('X', X).replaceAll('Y', Y).replaceAll('N', N).replaceAll('F', this.factorization.join(' * '));
};

function modInverse(a, m) {
  a = BigInt(a);
  m = BigInt(m);
  if (a <= 0n || a >= m || m <= 0n) {
    throw new RangeError();
  }
  // We use the extended Euclidean algorithm:
  let b = m;
  let [A, C] = [1n, 0n];
  while (b > 0n) {
    const q = a / b; // floor(a / b)
    [a, b] = [b, a - q * b];
    [A, C] = [C, A - q * C];
  }
  if (a !== 1n) {
    return 0n;
  }
  return A < 0n ? A + m : A;
}

function congruencesUsingContinuedFraction(primes, n) {
  const USE_LP_STRATEGY = true; // large primes
  let largePrimes = USE_LP_STRATEGY ? Object.create(null) : null; // prime -> congruence which needs this prime in base additionaly

  const lpStrategy = function (p, X, Y) {
    // https://ru.wikipedia.org/wiki/????????????????_??????????????#??????????????????_LP
    const lp = largePrimes[p];
    if (lp == undefined) {
      largePrimes[p] = {X: X, Y: Y};
    } else {
      const s = BigInt(p);
      const sInverse = modInverse(s, n);
      if (sInverse === 0n) {
        return new CongruenceOfsquareOfXminusYmoduloN(s, 0n, n, null);//?
      } else {
        const X1 = (sInverse * lp.X * X) % n;
        if (Y % s === 0n && lp.Y % s === 0n) {
          const Y1 = (lp.Y / s) * (Y / s);
          const factorization = getSmoothFactorization(Y1, primes);
          if (factorization != null) {
            return new CongruenceOfsquareOfXminusYmoduloN(X1, Y1, n, factorization);
          }
        }
      }
    }
    return null;
  };

  const primesProduct = product(primes);
  //const product1 = BigInt(primes.reduce((p, prime) => p * Number(prime) <= Number.MAX_SAFE_INTEGER ? p * Number(prime) : p, 1));
  const product1 = product(primes.slice(0, 11));
  n = BigInt(n);
  const d = ((n - n % 2n) / 2n);
  let [hprev, h] = [0n, 1n]; // previout and current convergent numerator modulo n
  const fraction = continuedFractionForSqrt(n);
  const iterator = {
    next: function congruencesUsingContinuedFraction() {
      let a = 0n;
      while ((a = fraction.next().value) != undefined) { // TODO: why do we stop after the first cycle ?
        // https://en.wikipedia.org/wiki/Continued_fraction#:~:text=The%20successive%20convergents%20are%20given%20by%20the%20formula
        [hprev, h] = [h, BigInt(a) * h + hprev];
        // optimization from the https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html :
        h = h % n;
        const X = h % n; // A_k mod n
        const Y = (X * X + d) % n - d; // (A_k)^2 mod n
        //console.log(X, Y);
        if (Y === 0n) {
          return {value: new CongruenceOfsquareOfXminusYmoduloN(X, 0n, n, null), done: false};
        } else {
          const s = isSmoothOverProduct(Y, primesProduct, product1);
          if (s === 1n) {
            return {value: new CongruenceOfsquareOfXminusYmoduloN(X, Y, n, getSmoothFactorization(Y, primes)), done: false};
          } else {
            if (USE_LP_STRATEGY) {
              const B = primes.length === 0 ? 1 : Number(primes[primes.length - 1]);
              const lp = Number(s);
              if (lp <= Math.min(B * B, Number.MAX_SAFE_INTEGER)) {
                // s is prime
                //if (!isPrime(s)) {
                //  throw new RangeError();
                //}
                const c = lpStrategy(lp, X, Y);
                if (c != null) {
                  return {value: c, done: false};
                }
              }
            }
          }
        }
      }
      return {value: undefined, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
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
BitSet.prototype.toggle = function (index) {
  if (index >= this.size) {
    throw new RangeError();
  }
  this.data[Math.floor(index / 32)] ^= (1 << (index % 32));
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
  return this.data.map(x => (x >>> 0).toString(2).padStart(32, '0').split('').reverse().join('')).join('').slice(0, this.size);
};

// pass factorizations with associated values to the next call
// returns linear combinations of vectors which result in zero vector by modulo 2
// (basis of the kernel of the matrix)
function solve(matrixSize) {
  // We build the augmented matrix in row-echelon form with permuted rows, which can grow up to matrixSize rows:
  const M = new Array(matrixSize).fill(null); // We will fill the matrix so pivot elements will be placed on the diagonal
  const associatedValues = new Array(matrixSize).fill(undefined);
  let nextSolution = null;
  let state = 1;
  const iterator = {
    next: function solve(tmp) {
      while (true) {
        if (state === 1) {
          state = 0;
          return {value: nextSolution, done: false};
        }
        state = 1;
        const [rawRow, associatedValue] = tmp;
        let row = new BitSet(matrixSize + matrixSize);
        for (let j = 0; j < rawRow.length; j += 1) {
          row.toggle(rawRow[j]);
        }
        // add row to the matrix maintaining it to be in row-echelon form:
        for (let pivotColumn = 0; pivotColumn < matrixSize && row != null; pivotColumn += 1) {
          if (row.has(pivotColumn)) {
            const pivotRow = M[pivotColumn];
            if (pivotRow != null) {
              // row-reduction:
              row.xor(pivotRow);
            } else {
              row.toggle(matrixSize + pivotColumn);
              associatedValues[pivotColumn] = associatedValue;
              M[pivotColumn] = row;
              row = null;
            }
          }
        }
        if (row != null) {
          // row has a solution
          // extract solution from the augmented part of the matrix:
          const solution = [];
          for (let i = 0; i < M.length; i += 1) {
            if (row.has(matrixSize + i)) {
              solution.push(associatedValues[i]);
            }
          }
          solution.push(associatedValue);
          nextSolution = solution;
        } else {
          nextSolution = null;
        }
      }
      //console.log(M.filter(x => x != null).map(x => x.toString()).join('\n'))
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
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

function abs(x) {
  return x < 0n ? -x : x;
}

function ContinuedFractionFactorization(N) {
  N = BigInt(N);
  //if (isPrime(N)) {
  //  return N;
  //}
  for (let k = 1n;; k += 1n) {
    const kN = k * N;
    // https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=optimal%20value :
    const B = Math.min(Math.floor(Math.sqrt(L(kN))), (1 << 25) - 1);
    // as we are check for smoothness the number P_k^2 - N * Q_k^2, where P_k/Q_k is a convergent,
    // if p is a factor of this number then P_k^2 - N * Q_k^2 = 0 (mod p),
    // after multiplication by Q_k^-2 we have (Q_k^-1*P_k)^2 = N (mod p),
    // x^2 = N (mod p), so N is a quadratic residue modulo N
    // ????. "??????????????????-???????????????? ???????????? ?? ????????????????????????" (??.??. ????????????, ??.??. ????????????????????), ???????????????? 202
    const primeBase = primes(B).filter(p => isQuadraticResidueModuloPrime(kN, p)).map(p => BigInt(p));
    const congruences = congruencesUsingContinuedFraction(primeBase, kN); // congruences X_k^2 = Y_k mod N, where Y_k is smooth over the prime base
    const solutions = solve(1 + primeBase.length); // find products of Y_k = Y, so that Y is a perfect square
    solutions.next();
    let c = null;
    while ((c = congruences.next().value) != undefined) {
      const solution = c.Y === 0n ? [c] : solutions.next([c.factorization.map(p => p === -1 ? 0 : 1 + primeBase.indexOf(p)), c]).value;
      if (solution != null) {
        const X = product(solution.map(c => c.X));
        const Y = product(solution.map(c => c.Y)); // = sqrt(X**2 % N)
        const x = X;
        const y = BigInt(sqrt(Y));
        console.assert(y * y === BigInt(Y));
        const g = gcd(abs(x + y), N);
        if (g !== 1n && g !== N) {
          return g;
        }
      }
    }
  }
}

ContinuedFractionFactorization.testables = {
  continuedFractionForSqrt: continuedFractionForSqrt,
  primes: primes,
  L: L,
  getSmoothFactorization: getSmoothFactorization,
  CongruenceOfsquareOfXminusYmoduloN: CongruenceOfsquareOfXminusYmoduloN,
  solve: solve,
  product: product,
  isQuadraticResidueModuloPrime: isQuadraticResidueModuloPrime
};

export default ContinuedFractionFactorization;
