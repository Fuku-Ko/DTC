program DTC
!Dicision Tree Classifier

implicit none
!class(must be integer)
integer, parameter :: cat_min = 0
integer, parameter :: cat_max = 2
!data number
integer, parameter :: nmax = 150
!raw number(last raw must be class)
integer, parameter :: mmax = 5
!exit criteria(data number by leaf)
integer, parameter :: min = 1
!max layer number
integer, parameter :: layermax = 10
integer, allocatable :: cat(:)
integer, allocatable :: y(:)
real, allocatable :: z(:,:,:)
real, allocatable :: x(:,:)
real, allocatable :: x_1(:,:)
integer, allocatable :: b(:)
real, allocatable :: c(:)
integer, allocatable :: d(:)
integer, allocatable :: e(:)
integer rightnum, leftnum, rightid, leftid, rightcat, leftcat
real right, both, left
integer i, j, k, l, lmax, maxnumber, sign, sign_cat
integer m, n, p, q, r, s, pmax, layer
real buf, score, threshold, scoremax
lmax = 2**layermax
allocate (z(1:lmax,1:nmax,1:mmax))
allocate (cat(cat_min:cat_max))
allocate (x(1:nmax,1:mmax))
allocate (b(1:lmax))
allocate (c(1:lmax))
allocate (d(1:lmax))
allocate (e(1:lmax))

open(80,file = 'iris-dataset.csv')
read(80, '()')
do n = 1, nmax
	read(80,*) x(n,:)
end do
close(80)

open(120,file = 'leaf.txt')
write(120,*)"layer leaf N data"

!initialize
b(1) = 0
c(1) = 0
d(1) = 0
e(1) = nmax
do i = 1, nmax
	do j = 1, mmax
		z(1,i,j) = x(i,j)
	end do
end do
do l = 2, lmax
	b(l) = 0
	c(l) = 0
	d(l) = 0
	e(l) = 0
	do i = 1, nmax
		do j = 1, mmax
			z(l,i,j) =0.0
		end do
	end do
end do

!by leaf (l)
do layer = 1, layermax
	do l = 2**(layer -1), 2**layer -1
		leftcat = 0
		rightcat = 0
		pmax = e(l)
		sign = 0
		!leaf categorize
		if (layer == layermax) then
			b(l) = 0
			c(l) = 0
		else
			!matrix prepare
			allocate (x_1(1:pmax,1:mmax))
			allocate (y(1:pmax))
			do i = 1, pmax
				do j = 1, mmax
					x_1(i,j) = z(l,i,j)
				end do
			end do
			!by variable (search best threshold)
			leftid = 0
			rightid = 0
			do m = 1, mmax-1
				!sort
				do i = 1, pmax-1
					do j = i+1, pmax
						if (x_1(i,m)>x_1(j,m)) then
							do k = 1, mmax
								buf = x_1(i,k)
								x_1(i,k) = x_1(j,k)
								x_1(j,k) = buf
							end do
						end if
					end do
				end do
				do i = 1, pmax
					y(i) = nint(x_1(i, mmax))
				end do
				!score calculate and threshold dicision
				if(pmax>1)then
					do i = 1, pmax-1
						if (x_1(i,m) == x_1(i+1,m)) then
						else
							left = 0.0
							right = 0.0
							both = 0.0
							score = 0.0
							do j = cat_min, cat_max
								cat(j) = 0
							end do
							do j = 1, i
								s = y(j)
								cat(s) = cat(s) +1
							end do
							do j = cat_min, cat_max
								if(cat(j) == 0)then
								else
									left = left -(cat(j)*1.0/(i*1.0))*log(cat(j)*1.0/(i*1.0))
								end if
							end do
							leftnum = 0
							leftid = 0
							do j = cat_min, cat_max
								if(leftnum<cat(j))then
									leftnum = cat(j)
									leftid = j
								end if
							end do
							do j = cat_min, cat_max
								cat(j) = 0
							end do
							do j = i +1, pmax
								s = y(j)
								cat(s) = cat(s) +1
							end do
							do j = cat_min, cat_max
								if(cat(j) == 0)then
								else
									right = right -(cat(j)*1.0/((pmax -i)*1.0))*log(cat(j)*1.0/((pmax -i)*1.0))
								end if
							end do
							rightnum = 0
							rightid = 0
							do j = cat_min, cat_max
								if(rightnum<cat(j))then
									rightnum = cat(j)
									rightid = j
								end if
							end do
							do j = cat_min, cat_max
								cat(j) = 0
							end do
							do j = 1, pmax
								s = y(j)
								cat(s) = cat(s) +1
							end do
							do j = cat_min, cat_max
								if(cat(j) == 0)then
								else
									both = both -(cat(j)*1.0/(pmax*1.0))*log(cat(j)*1.0/(pmax*1.0))
								end if
							end do
							score = both -left*i/(pmax*1.0) -right*(pmax-i)/(pmax*1.0)
							if(score>scoremax .or. sign == 0)then
								sign = 1
								scoremax = score
								maxnumber = m
								leftcat = leftid
								rightcat = rightid
								threshold = (x_1(i,m) + x_1(i+1,m))*0.5
							end if
						end if
					end do
				end if
			end do
			
			!not make leaf(all same value or category)
			sign_cat = 0
			do i = cat_min, cat_max
				cat(i) = 0
			end do
			do i = 1, pmax
				s = y(i)
				cat(s) = cat(s) +1
			end do
			do i = cat_min, cat_max
				if(cat(i) == 0 .or. cat(i) == pmax)then
				else
					sign_cat = 1
				end if
			end do
			if (sign == 0 .or. sign_cat == 0) then
				if (sign_cat == 0) then
					do i = 1, pmax
						write(120,*) layer, l, e(l), x_1(i,:)
					end do
				end if
				b(l) = 0
				c(l) = 0
			else
				!make leaf by threshold
				q = 1
				r = 1
				do i = 1, pmax
					if (x_1(i,maxnumber)<threshold) then
						do j = 1, mmax
							z(2*l,q,j) = x_1(i,j)
						end do
						e(2*l) = e(2*l) +1
						q = q +1
					else
						do j = 1, mmax
							z(2*l +1,r,j) = x_1(i,j)
						end do
						e(2*l + 1) = e(2*l + 1) +1
						r = r +1
					end if
				end do
				
				!memorize leaf
				if(e(2*l)<min .or. e(2*l + 1)<min)then
					b(l) = 0
					c(l) = 0
					e(2*l) = 0
					e(2*l +1) = 0
				else 
					b(l) = maxnumber
					c(l) = threshold
					d(2*l) = leftcat
					d(2*l + 1) = rightcat
				end if
			end if
			deallocate (x_1)
			deallocate (y)
		end if
	end do
end do

close(120)

open(100,file = 'model.txt')
open(110,file = 'check.txt')
write(100,*) 'leaf, ','category(0 means stop), ','threshold, ','class(predict), ','data_number'
do layer = 1, layermax
	write(110,*)"layer=", layer
	write(110,*) 'leaf, ','category(0 means stop), ','threshold, ','class(predict), ','data_number'
	do l = 2**(layer -1), 2**layer -1
		write(100,*) l, b(l), c(l), d(l), e(l)
		if(e(l)>0)then
			write(110,*) l, b(l), c(l), d(l), e(l)
		end if
	end do
end do

close(100)
close(110)

deallocate (z)
deallocate (cat)
deallocate (x)
deallocate (b)
deallocate (c)
deallocate (d)
deallocate (e)

end program DTC
