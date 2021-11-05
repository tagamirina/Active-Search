#!/bin/bash

for num in {1..110}

do
	for tmp in {1..10}
	do

	if [ "$num" -ge 100 ]

	then

		if [ "$tmp" -ge 10 ]

		then
			./active.exe -i image/tmp0$tmp.png -ii image/img$num.png -map actmap/map0${tmp}_$num.png -ss result/result0${tmp}_$num.txt -res result/matchingresults0${tmp}_$num.png
		else
			./active.exe -i image/tmp00$tmp.png -ii image/img$num.png -map actmap/map00${tmp}_$num.png -ss result/result00${tmp}_$num.txt -res result/matchingresults00${tmp}_$num.png              
		fi

	else

		if [ "$num" -ge 10 ]

		then

			if [ "$tmp" -ge 10 ]

			then
				./active.exe -i image/tmp0$tmp.png -ii image/img0$num.png -map actmap/map0${tmp}_0$num.png -ss result/result0${tmp}_0$num.txt -res result/matchingresults0${tmp}_0$num.png
			else
				./active.exe -i image/tmp00$tmp.png -ii image/img0$num.png -map actmap/map00${tmp}_0$num.png -ss result/result00${tmp}_0$num.txt -res result/matchingresults00${tmp}_0$num.png
			fi

		else

			if [ "$tmp" -ge 10 ]

			then
				./active.exe -i image/tmp0$tmp.png -ii image/img00$num.png -map actmap/map0${tmp}_00$num.png -ss result/result0${tmp}_00$num.txt -res result/matchingresults0${tmp}_00$num.png
			else
				./active.exe -i image/tmp00$tmp.png -ii image/img00$num.png -map actmap/map00${tmp}_00$num.png -ss result/result00${tmp}_00$num.txt -res result/matchingresults00${tmp}_00$num.png
			fi

		fi

	fi

	done
done