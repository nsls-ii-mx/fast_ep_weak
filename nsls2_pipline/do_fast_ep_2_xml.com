echo '<root>' > fast_ep.xml
echo '<OutputData>' >> fast_ep.xml

echo '<BestHand>' >> fast_ep.xml
grep Best fast_ep.log | grep hand | sed 's/Best hand://g' >> fast_ep.xml
echo '</BestHand>' >> fast_ep.xml

echo '<NSites>' >> fast_ep.xml
grep Best fast_ep.log | grep nsites | sed 's/Best nsites://g' >> fast_ep.xml
echo '</NSites>' >> fast_ep.xml

echo '</OutputData>' >> fast_ep.xml
echo '</root>' >> fast_ep.xml
