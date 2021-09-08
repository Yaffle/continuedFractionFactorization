# continuedFractionFactorization
Continued Fraction Factorization (https://en.wikipedia.org/wiki/Continued_fraction_factorization) in JavaScript using native BigInt

There is description at https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html . See also links in the code.

Example
```javascript
import factorize from './continuedFractionFactorization.js';
console.time();
const f = factorize(2n**128n + 1n);
console.timeEnd();
// ~13 seconds
console.assert(f === 5704689200685129054721n || f === 59649589127497217n, f);
```
