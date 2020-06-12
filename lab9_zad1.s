
	.data
	
control_reg: 	.word	0			# Control Register
number:			.tfloat	2.0			# x - 80 bitow
napis:			.string "Value = %1.25Lf\n"	
	.equ precision,		0xFCFF 
	.equ single,		0x0000
	.equ DOUBLE,		0x0200
	.equ LONG_DOUBLE,	0x0300

	.equ round,				0xF3FF
	.equ nearest_even,		0x0000
	.equ round_down,		0x0400
	.equ round_up,			0x0800
	.equ truncate,			0x0C00


	.equ flags,				0xFFE0

	.equ zerodivide,			0x003B
	.equ invalid_oper,			0x003E
	

	.text
	.global main

	
main:
	sub		$8,%rsp    # przesuniecie ramki stosu

	finit			# inicjalizacja FPU,  rejestry są czyszczone
	
	fstcw control_reg	# zawartosc rejestru control register zostaje wczytana z ramu do zmiennej

	mov control_reg, %ax	# control register do rejstru %ax - będziemy go manipulować flagami

	and $precision, %ax	# ustawiamy precyzję, bity PC w control register (odpowiadające za precyzję) są zerowane
#	or $single, %ax		# stan bitóW PC - 00
#	or $DOUBLE, %ax		# stan bitóW PC - 10
	or $LONG_DOUBLE, %ax # stan bitóW PC - 11

	and $round, %ax	# ustawienie zaokrąglenia
#	or $nearest_even, %ax
#	or $round_down, %ax
#	or $round_up, %ax
	or $truncate, %ax

	and $flags,%ax # ustawienie wyjątków
	or $zerodivide, %ax
#	or $invalid_oper, %ax


	mov %ax, control_reg	# %ax to variable

	fldcw control_reg	# ładuje zawartość zmiennej do controlregister  

zerodivide:
	fldz
	fld1
	fdivp

denormal_operation:
#	fld1
#	fchs
#	fsqrt

invalid:
#	fadd	%st(1),%st(0)



stackoverflow:
#	fldpi
#	fldpi
#	fldpi
#	fldpi
#	fldpi
#	fldpi
#	fldpi
#	fldpi



	fldt number			
	fsqrt			# sqrt ST(0)  
	fldt number			# x -> ST(0)
	fsqrt			# sqrt(x) -> ST(0) 
	fmulp			# sqrt(x)^2 -> ST(0) | mulitply and pop
	fldt number			# x -> ST(0)	x,sqrt(x)*sqrt(x)
	fsubrp			# sqrt(x)^2 - x

	sub $16, %rsp		# Potrzeba 10 bajtow miejsca trzeba zapewnic wyrownanie stosu - 8 bajtow za malo, wiec robi sie miejsce dla 16 bajtow
	fstpt (%rsp)		# store in rsp and pop from stack
	mov $napis, %rdi	# zgodnie z abi
	xor %al, %al		# nie są używane rejestry wektorowe
	call printf	
	add $16, %rsp		# przesuwam stos spowrotem


	add 	$8,%rsp		# wyrownuje ramke stosu
	ret












