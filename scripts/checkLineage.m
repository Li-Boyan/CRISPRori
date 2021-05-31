close all
lineage.timeCourse(lineage.timeCourse.Time==0,:) = [];
repinit_idx = find(lineage.timeCourse.OriCNum(1:(end-1)) < lineage.timeCourse.OriCNum(2:end));
figure, yyaxis left, hold on
plot(lineage.timeCourse.Time, lineage.timeCourse.Total_proteins/3e6)
plot(lineage.timeCourse.Time(repinit_idx), lineage.timeCourse.Total_proteins(repinit_idx)/3e6, 'r.')
yyaxis right
plot(lineage.timeCourse.Time, lineage.timeCourse.GenomeNum)
figure, hold on
plot(lineage.timeCourse.Time, lineage.timeCourse.DnaA)
plot(lineage.timeCourse.Time, lineage.timeCourse.DnaAatp)
plot(lineage.timeCourse.Time, lineage.timeCourse.DnaAatp_free)
plot(lineage.timeCourse.Time, lineage.timeCourse.DnaAadp)
legend({'DnaA', 'DnaA-ATP', 'DnaA-ATP Free', 'DnaA-ADP'})
hold off