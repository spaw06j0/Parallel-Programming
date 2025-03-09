#include "PPintrin.h"

// implementation of absSerial(), but it is vectorized using PP intrinsics
void absVector(float *values, float *output, int N)
{
  __pp_vec_float x;
  __pp_vec_float result;
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {

    // All ones
    maskAll = _pp_init_ones();

    // All zeros
    maskIsNegative = _pp_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _pp_vload_float(x, values + i, maskAll); // x = values[i];

    // Set mask according to predicate
    _pp_vlt_float(maskIsNegative, x, zero, maskAll); // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _pp_vsub_float(result, zero, x, maskIsNegative); //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _pp_mask_not(maskIsNegative); // } else {

    // Execute instruction ("else" clause)
    _pp_vload_float(result, values + i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _pp_vstore_float(output + i, result, maskAll);
  }
}

void clampedExpVector(float *values, int *exponents, float *output, int N)
{
  //
  // PP STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    __pp_mask maskAll, maskIsZero, maskIsNotZero;
    int remainElements = N - i;
    if (remainElements < VECTOR_WIDTH) {
      maskAll = _pp_init_ones(remainElements);
    } else {
      maskAll = _pp_init_ones();
    }
    __pp_vec_float x;
    __pp_vec_int y;
    __pp_vec_int count;
    __pp_vec_float result;

    __pp_vec_int zeros = _pp_vset_int(0);
    __pp_vec_int ones = _pp_vset_int(1);
    __pp_vec_float max = _pp_vset_float(9.999999f);
    // load values and exponents
    _pp_vload_float(x, values + i, maskAll);
    _pp_vload_int(y, exponents + i, maskAll);
    // if y == 0 result = 1.f
    _pp_veq_int(maskIsZero, y, zeros, maskAll);
    _pp_vset_float(result, 1.f, maskIsZero);

    // Handle non-zero exponents
    _pp_vgt_int(maskIsNotZero, y, zeros, maskAll);
    _pp_vmove_float(result, x, maskIsNotZero);  // Initialize result = x for non-zero exponents
    _pp_vsub_int(count, y, ones, maskIsNotZero); // count = y - 1
    _pp_vgt_int(maskIsNotZero, count, zeros, maskAll);
    // Compute x^y for non-zero exponents
    while (_pp_cntbits(maskIsNotZero) > 0) {
        _pp_vmult_float(result, result, x, maskIsNotZero); // result *= x
        _pp_vsub_int(count, count, ones, maskIsNotZero);   // count--
        _pp_vgt_int(maskIsNotZero, count, zeros, maskAll); // Update mask for count > 0
    }

    // Clamp results to 9.999999f
    __pp_mask maskClamp;
    _pp_vgt_float(maskClamp, result, max, maskAll);
    _pp_vset_float(result, 9.999999f, maskClamp);

    _pp_vstore_float(output + i, result, maskAll);
  }
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N)
{
  __pp_vec_float x;
  __pp_vec_float sum = _pp_vset_float(0.f);
  __pp_mask maskAll = _pp_init_ones();
  
  // Accumulate all chunks
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    _pp_vload_float(x, values + i, maskAll);
    _pp_vadd_float(sum, sum, x, maskAll);
  }
  // [0 1 2 3 4 5 6 7] -> [0 2 4 6 1 3 5 7]
  // [0+1, 0+1, 2+3, 2+3, 4+5, 4+5, 6+7, 6+7]
  // [0+1, 2+3, 4+5, 6+7, 0+1, 2+3, 4+5, 6+7]
  // [0+1+2+3, 0+1+2+3, 4+5+6+7, 4+5+6+7]
  // size 2 -> 1, size 4 -> 2, size 8 -> 3
  // Reduction phase using the correct pattern
  int shiftCounter = VECTOR_WIDTH;
  while(shiftCounter != 1) {
    _pp_hadd_float(sum, sum);
    _pp_interleave_float(sum, sum);
    shiftCounter >>= 1;
  }
  
  // Return the first element which now contains the sum
  return sum.value[0];
}