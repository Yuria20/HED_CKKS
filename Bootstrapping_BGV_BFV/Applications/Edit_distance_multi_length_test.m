load "Applications/Edit_distance.m";

// Basic checks
assert scheme eq "BFV";
assert GetLTVersion() eq 4;
assert (r eq 1) and (e eq 2);
assert IsPowerOfTwo(n_prime) and IsPowerOfTwo(n_double_prime);
assert (n_prime le n) and (n_double_prime lt n_prime);
assert &*gbfv_mat_dimensions eq Minimum(n_prime, n_double_prime ^ 2);
assert (#gbfv_mat_dimensions eq 0) or (Minimum(gbfv_mat_dimensions) gt 1);

// Check if the used moduli are in the cyclotomic ring
assert &and[ToCyclotomicField(allModuli[index + 1]) * InvertOverField(allModuli[index]) in MaximalOrder(cyclo_field) :
            index in [1..#allModuli - 1]];

// Check if decomposition parameters are consistent with intermediate moduli
assert #intModuli lt Maximum(#gbfv_mat_dimensions, 1);
assert gbfvExponent * n_prime eq n_double_prime * n;
assert (#gbfv_mat_dimensions eq 0) or
       (&and[intExponents[index] * &*gbfv_mat_dimensions[#gbfv_mat_dimensions - #intExponents + index..#gbfv_mat_dimensions] *
             Maximum(n_prime div (n_double_prime ^ 2), 1) ge n : index in [1..#intExponents]]);

// Only full packing case is implemented
assert n eq n_prime;



// Plaintext edit distance (Levenshtein) for verification
function PlainEditDistance(a, b)
    prev := [0..#b];
    for i := 1 to #a do
        cur := [i];
        for j := 1 to #b do
            Append(~cur, Minimum([prev[j + 1] + 1, cur[j] + 1, prev[j] + ((a[i] eq b[j]) select 0 else 1)]));
        end for;
        prev := cur;
    end for;
    return prev[#b + 1];
end function;

// Experiment settings
// Each entry is one string pair; both strings in a pair must have the same length
// (powers of two give clean timings)
test_pairs := [
    <"KITTENSS", "SITTINGS">,                                               // length 8
    <"INTENTIONEXECUTE", "EXECUTIONINTENTS">,                               // length 16
    <"THEQUICKBROWNFOXJUMPSOVERLAZYDOG", "THEQUIETBROWNCATJUMPEDOVERTHEDOG"> // length 32
];
alphabet_size_encrypt := 36;
alphabet_size_evaluate := 128;

// Length-independent key generation (done only once)
sk, pk := GenKeyPair(); rk := GenRelinKey(sk);
recrypt_variables := GenerateGBFVBetterRecryptVariables(sk, pk, B);



for pair in test_pairs do
    str_a := [pair[1]];
    str_b := [pair[2]];
    len := #pair[1];
    assert #pair[2] eq len;
    nb_pairs := #str_a;
    label := "len " cat IntegerToString(len);

    // Length-dependent rotation keys
    len_str := CeilPowerOfTwo(len);
    expand_keys, extract_keys, key_right, key_left, boot_shift_key := PrecompKeys(sk, len_str);

    // Encryption, edit distance computation and decryption
    ctxt_a, ctxt_b := EncryptStrings(str_a, str_b, alphabet_size_encrypt, pk);
    timer := StartTiming();
    ctxt_edit := EditDistance(ctxt_a, ctxt_b, len_str, alphabet_size_evaluate, recrypt_variables, expand_keys, extract_keys,
                              key_right, key_left, boot_shift_key, rk: iterations_per_bootstrap := 2);
    StopTiming(timer: message := "edit distance " cat label);

    // Compare with plaintext edit distance
    result := DecryptEditDistances(ctxt_edit, sk, nb_pairs);
    expected := [Z | PlainEditDistance(str_a[i], str_b[i]) : i in [1..nb_pairs]];
    print label, "homomorphic:", result, "plaintext:", expected, "match:", result eq expected;

    // Check correctness of SEAL output
    expected_result := ScaleAndRound(Decrypt(ctxt_edit, sk), p, gbfvModulus) mod p;
    _ := Decrypt(SetPlaintextModulus(ctxt_edit, p), sk: print_result := false, check_correctness := true,
                                                        expected_result := expected_result, message := label);
end for;

load "Traces/Post_processing.m";
