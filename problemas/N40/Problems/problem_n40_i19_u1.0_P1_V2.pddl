(define (problem compiledcode)
(:domain quantum-gate)


(:objects
   q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21 q22 q23 q24 q25 q26 q27 q28 q29 q30 q31 q32 q33 q34 q35 q36 q37 q38 q39 q40  - qstate)


(:init
   (located_at_1 q1)
   (located_at_2 q2)
   (located_at_3 q3)
   (located_at_4 q4)
   (located_at_5 q5)
   (located_at_6 q6)
   (located_at_7 q7)
   (located_at_8 q8)
   (located_at_9 q9)
   (located_at_10 q10)
   (located_at_11 q11)
   (located_at_12 q12)
   (located_at_13 q13)
   (located_at_14 q14)
   (located_at_15 q15)
   (located_at_16 q16)
   (located_at_17 q17)
   (located_at_18 q18)
   (located_at_19 q19)
   (located_at_20 q20)
   (located_at_21 q21)
   (located_at_22 q22)
   (located_at_23 q23)
   (located_at_24 q24)
   (located_at_25 q25)
   (located_at_26 q26)
   (located_at_27 q27)
   (located_at_28 q28)
   (located_at_29 q29)
   (located_at_30 q30)
   (located_at_31 q31)
   (located_at_32 q32)
   (located_at_33 q33)
   (located_at_34 q34)
   (located_at_35 q35)
   (located_at_36 q36)
   (located_at_37 q37)
   (located_at_38 q38)
   (located_at_39 q39)
   (located_at_40 q40)
)


(:goal
   (and
      (U_GOAL q22 q33)
      (U_GOAL q11 q31)
      (U_GOAL q20 q30)
      (U_GOAL q2 q9)
      (U_GOAL q12 q25)
      (U_GOAL q12 q36)
      (U_GOAL q4 q10)
      (U_GOAL q27 q28)
      (U_GOAL q16 q32)
      (U_GOAL q24 q30)
      (U_GOAL q11 q13)
      (U_GOAL q24 q35)
      (U_GOAL q25 q39)
      (U_GOAL q7 q34)
      (U_GOAL q9 q28)
      (U_GOAL q33 q40)
      (U_GOAL q6 q23)
      (U_GOAL q11 q23)
      (U_GOAL q21 q23)
      (U_GOAL q4 q27)
      (U_GOAL q18 q40)
      (U_GOAL q3 q19)
      (U_GOAL q20 q35)
      (U_GOAL q8 q36)
      (U_GOAL q11 q25)
      (U_GOAL q17 q38)
      (U_GOAL q20 q39)
      (U_GOAL q12 q26)
      (U_GOAL q10 q31)
      (U_GOAL q14 q34)
      (U_GOAL q25 q26)
      (U_GOAL q7 q15)
      (U_GOAL q4 q14)
      (U_GOAL q20 q21)
      (U_GOAL q27 q30)
      (U_GOAL q32 q38)
      (U_GOAL q7 q26)
      (U_GOAL q9 q16)
      (U_GOAL q3 q9)
      (U_GOAL q22 q38)
   )
)


(:metric minimize (total-time))
)
