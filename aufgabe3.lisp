(defun my-xor (n1 n2)
  (if (= n1 n2) 0 1))


(defun binary-to-graycode (binary)
  (loop with size = (length binary)
        with graycode = (make-array size :element-type 'bit :initial-element (bit binary (1- size)))
          for i from (- size 2) downto 0
          do (setf (bit graycode i) (my-xor (bit binary (1+ i)) (bit binary i)))
        finally (return graycode)))


(defun graycode-to-binary (graycode)
  (loop with size = (length graycode)
        with binary = (make-array size :element-type 'bit :initial-element (bit graycode (1- size)))
          for i from (- size 2) downto 0
          do (setf (bit binary i) (my-xor (bit binary (1+ i)) (bit graycode i)))
        finally (return binary)))


(defun random-bits (n)
  (loop with bit-arr = (make-array n :element-type 'bit)
          for i below n 
          do (setf (bit bit-arr i) (random 2))
        finally (return bit-arr)))


(defun init-population (amount gene-length)
  (loop with population = (make-array amount :element-type 'bit-vector)
          for i below amount and 
          for genes = (random-bits gene-length)
            do (setf (aref population i) genes)
          finally (return population)))
  

(defun mutate (bits-to-mutate)
  (let ((bit-mask (make-array (length bits-to-mutate) :element-type 'bit)))
    (setf (bit bit-mask (random (length bits-to-mutate))) 1)
    (bit-xor bits-to-mutate bit-mask)))


;; sets all bits to zero, only the bits between start and opt-end keep its state
;;  if opt-end is omitted, bits from start till (length bits) will be extracted
(defun extract-bit-sequence (bits start &optional opt-end)
  (let ((bit-mask (make-array (length bits) :element-type 'bit))
        (end (if opt-end opt-end (length bits))))
      (loop for i from start below end 
              do (setf (bit bit-mask i) 1))
      (bit-and bit-mask bits)))


;; indices must be in ascending order
(defun crossover (bit-mother bit-father split-indices)
;; Helper function
  (labels ((crossover-tr (bit-mom bit-dad index-lst acc)
              (if (null index-lst)
                  acc
                  (crossover-tr bit-dad 
                                bit-mom 
                                (rest index-lst)
                                (bit-ior (extract-bit-sequence bit-mom 
                                                               (first index-lst) 
                                                               (second index-lst)) 
                                          acc)))))
    (values (crossover-tr bit-mother 
                          bit-father
                          (if (= 0 (first split-indices)) split-indices (cons 0 split-indices))
                          (make-array (length bit-mother) :element-type 'bit))
            (crossover-tr bit-father 
                          bit-mother
                          (if (= 0 (first split-indices)) split-indices (cons 0 split-indices))
                          (make-array (length bit-mother) :element-type 'bit)))))


(defun random-crossover-indices (ct max)
  (sort (loop for i below ct collect (random max)) #'<))


(defun select-hypo (hypos prob-function)
  (loop with size  = (length hypos) and 
             sum   = 0 and
             rand  = (random 1.0)
        with index = (random size)
          do 
            (setf index (mod (1+ index) size))
            (setf sum (+ sum (funcall prob-function (aref hypos index))))
          when (< rand sum) return index))


(defun selection (hypos prob-function ct dest &key (offset 0))
  (loop for i below ct 
          do (setf (aref dest (+ i offset)) (aref hypos (select-hypo hypos prob-function)))))


(defun fitness (gene volumes constant)
  (flet ((square (n) (* n n)))
    (handler-case (exp (* (- constant) 
                          (square (- 100 
                                  (loop for i below (length gene) 
                                          sum (* (bit gene i) (aref volumes i)))))))
                  (floating-point-underflow nil 0.0))))


(defun create-prob-function (population fitness-func)
  (let ((sum-of-fitness (reduce #'(lambda (acc gene) (+ acc (funcall fitness-func gene)))
                                population
                                :initial-value 0)))
    #'(lambda (gene) (/ (funcall fitness-func gene) sum-of-fitness))))


(defun find-max-hypo (population fitness-func)
  (loop with max = 0
          for i below (length population)
            when (< (funcall fitness-func (aref population max))
                    (funcall fitness-func (aref population i)))
                  do (setf max i)
          finally (return max)))


(defun calc-volume (gene volumes)
  (loop for x below (length gene) sum (* (aref gene x) (aref volumes x))))


(defun evolutionary (p r m gene-length run-ct crossover-indices-ct fitness-func)
  (loop with population    = (init-population p gene-length)
        and  next-gen      = (make-array p :element-type 'bit-vector)
        and  mutations     = (round (* m p))
        and  crossovers    = (round (* r p))

        for prob-function  = (create-prob-function population fitness-func) and
        for best-hypo      = (aref population (find-max-hypo population fitness-func)) and
        for i below run-ct
        do
          (format t "~%~A~%" (funcall fitness-func best-hypo))
          ;; (format t "~%~A~%" best-hypo)
          (selection population prob-function (- p crossovers) next-gen)
          (selection population prob-function crossovers next-gen :offset (- p crossovers))
          ;; crossover
          (loop for i from (- p crossovers) below (1- p) by 2
            do (setf (values (aref next-gen i) (aref next-gen (1+ i)))
                     (crossover (aref next-gen i) 
                                (aref next-gen (1+ i)) 
                                ;; '(20 40 60 80))))
                                (random-crossover-indices crossover-indices-ct gene-length))))
          ;; Mutate
          (loop for index = (random p)
                for i below mutations
                  do (setf (aref next-gen index) (mutate (aref next-gen index))))
          ;; Swap the array's roles
          (let ((tmp population))
            (setf population next-gen)
            (setf next-gen tmp))
        finally (return (values best-hypo next-gen))))


(defvar *volumes* (loop with x = (make-array 100) for i below 100 do (setf (aref x i) (1+ (random 10.0))) finally (return x)))

(defvar *constant* 0.0001)

(defvar *fitness-function* #'(lambda (gene) (fitness gene *volumes* *constant*)))
(defvar a)
(defvar b)
(setf (values a b) (evolutionary 10 0.4 0.3 100 50 1 *fitness-function*))

(format t "Resulting fitness:~%~A~%" (funcall *fitness-function* a))
;; (format t "Resulting gene:~%~A~%" a)
(format t "Resulting volume:~%~A~%" (loop for i below (length a) sum (* (aref a i) (aref *volumes* i))))
