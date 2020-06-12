
	.data
	
counter:      .quad	5000000

control_reg:	.word	0	
napis:			.string "Pi = %1.25Lf\n"	

    .equ precision,		0xFCFF 
	.equ single,		0x0000
	.equ DOUBLE,		0x0200
	.equ LONG_DOUBLE,	0x0300




	.text
	.global main

	
main:
	push %rbp

	finit			# x87 initialzation
	
	fstcw control_reg	# CR -> variable

	mov control_reg, %bx	# CR to %ax

	and $precision, %bx	# set precision
#	or $single, %ax
#	or $DOUBLE, %ax
	or $LONG_DOUBLE, %bx

	mov %bx, control_reg	

	fldcw control_reg	

	movq counter, %rcx	

# Przygotowanie liczb 
	fld1			# 1
	fld1			# 1, 1
	fadd %st(1), %st(0)	# 2, 1
	fld %st(0)		# 2, 2, 1
	fadd %st(1), %st(0)	# 4, 2, 1
	fchs			# -4, 2, 1
	fldz			# 0, -4, 2, 1	

next:
	fld %st(1)		# -4, 0, -4, 2, 1
	fchs			# 4, 0, -4, 2, 1
	fst %st(2)		# 4, 0, 4, 2, 1
	fdiv %st(4)		# 4/1, 0, 4, 2, 1
	faddp			# 4/1, 4, 2, 1
	fld %st(3)		# 1, 4/1, 4, 2, 1
	fadd %st(3), %st(0)	# 3, 4/1, 4, 2, 1
	fstp %st(4)		# 4/1, 4, 2, 3


	dec %rcx
	jnz next
	
	sub $16, %rsp		
	fstpt (%rsp)		
	mov $napis, %rdi	
	xor %al, %al		
	call printf	
	add $16, %rsp		


	pop %rbp

	ret
