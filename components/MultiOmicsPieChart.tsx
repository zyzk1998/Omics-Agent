/**
 * MultiOmicsPieChart - Global Multi-omics Market Share (2025 Forecast)
 *
 * Dependencies: npm install recharts
 * Styling: Tailwind CSS
 */
/**
 * MultiOmicsPieChart - Global Multi-omics Market Share (2025 Forecast)
 * Dependencies: npm install recharts | Styling: Tailwind CSS
 */
'use client';

import React, { useState, useEffect } from 'react';
import {
  PieChart,
  Pie,
  Cell,
  Tooltip,
  Legend,
  ResponsiveContainer,
  Sector,
} from 'recharts';

// Data with Key Drivers for tooltip
const MARKET_DATA = [
  {
    name: 'Genomics',
    value: 42,
    color: '#3B82F6',
    keyDriver: 'Clinical WGS & Oncology',
  },
  {
    name: 'Proteomics',
    value: 24,
    color: '#10B981',
    keyDriver: 'Biomarker Discovery & Therapeutics',
  },
  {
    name: 'Transcriptomics',
    value: 18,
    color: '#8B5CF6',
    keyDriver: 'Single-cell & Spatial Transcriptomics',
  },
  {
    name: 'Metabolomics',
    value: 10,
    color: '#F59E0B',
    keyDriver: 'Precision Nutrition & Drug Metabolism',
  },
  {
    name: 'Epigenomics',
    value: 4,
    color: '#EF4444',
    keyDriver: 'Gene Regulation & Disease Mechanisms',
  },
  {
    name: 'Microbiomics & Others',
    value: 2,
    color: '#6B7280',
    keyDriver: 'Gut Microbiome & Host-Microbe Interactions',
  },
];

// Custom active shape with expanded sector on hover
const renderActiveShape = (props: {
  cx: number;
  cy: number;
  midAngle: number;
  innerRadius: number;
  outerRadius: number;
  startAngle: number;
  endAngle: number;
  fill: string;
  payload: (typeof MARKET_DATA)[0];
  isActive: boolean;
}) => {
  const {
    cx,
    cy,
    innerRadius,
    outerRadius,
    startAngle,
    endAngle,
    fill,
    isActive,
  } = props;

  const expandRadius = isActive ? 8 : 0;

  return (
    <Sector
      cx={cx}
      cy={cy}
      innerRadius={innerRadius}
      outerRadius={outerRadius + expandRadius}
      startAngle={startAngle}
      endAngle={endAngle}
      fill={fill}
      style={{
        filter: isActive ? 'drop-shadow(0 4px 6px rgba(0,0,0,0.1))' : undefined,
        transition: 'all 0.3s ease-out',
      }}
    />
  );
};

// Custom Tooltip
const CustomTooltip = ({
  active,
  payload,
}: {
  active?: boolean;
  payload?: Array<{ payload: (typeof MARKET_DATA)[0] }>;
}) => {
  if (!active || !payload?.length) return null;

  const data = payload[0].payload;

  return (
    <div className="rounded-lg border border-slate-200 bg-white px-4 py-3 shadow-lg">
      <div className="space-y-1">
        <p className="font-semibold text-slate-800">{data.name}</p>
        <p className="text-sm text-slate-600">
          <span className="font-medium">{data.value}%</span> Market Share
        </p>
        <p className="text-xs text-slate-500">
          <span className="font-medium">Driver:</span> {data.keyDriver}
        </p>
      </div>
    </div>
  );
};

// Custom Legend
const renderLegend = (props: {
  payload?: Array<{ value?: string; payload?: { name: string }; color: string }>;
}) => {
  const { payload } = props;
  if (!payload) return null;

  return (
    <div className="flex flex-wrap justify-center gap-x-6 gap-y-2 pt-4">
      {payload.map((entry) => {
        const label = entry.value ?? entry.payload?.name ?? '';
        return (
          <div key={label} className="flex items-center gap-2">
            <span
              className="h-3 w-3 shrink-0 rounded-full"
              style={{ backgroundColor: entry.color }}
            />
            <span className="text-sm font-medium text-slate-600">{label}</span>
          </div>
        );
      })}
    </div>
  );
};

export default function MultiOmicsPieChart() {
  const [activeIndex, setActiveIndex] = useState<number | undefined>(undefined);
  const [mounted, setMounted] = useState(false);

  useEffect(() => {
    setMounted(true);
  }, []);

  return (
    <div
      className={`w-full max-w-2xl transition-all duration-500 ease-out ${
        mounted ? 'translate-y-0 opacity-100' : 'translate-y-4 opacity-0'
      }`}
    >
      <div className="rounded-xl border border-slate-200 bg-white p-6 shadow-sm">
        <h3 className="mb-6 text-center text-lg font-semibold text-slate-800">
          Global Multi-omics Market Share (2025 Forecast)
        </h3>

        <ResponsiveContainer width="100%" height={320}>
          <PieChart margin={{ top: 10, right: 10, bottom: 10, left: 10 }}>
            <Pie
              data={MARKET_DATA}
              cx="50%"
              cy="50%"
              innerRadius={60}
              outerRadius={100}
              paddingAngle={1}
              dataKey="value"
              nameKey="name"
              isAnimationActive
              animationBegin={200}
              animationDuration={800}
              onMouseEnter={(_, index) => setActiveIndex(index)}
              onMouseLeave={() => setActiveIndex(undefined)}
              activeIndex={activeIndex}
              activeShape={(props) =>
                renderActiveShape({
                  ...props,
                  isActive: true, // activeShape only renders for hovered sector
                })
              }
            >
              {MARKET_DATA.map((entry, index) => (
                <Cell
                  key={entry.name}
                  fill={entry.color}
                  stroke="white"
                  strokeWidth={2}
                  style={{
                    opacity: activeIndex === undefined || activeIndex === index ? 1 : 0.7,
                    transition: 'opacity 0.2s ease',
                  }}
                />
              ))}
            </Pie>

            <Tooltip content={<CustomTooltip />} />
            <Legend content={renderLegend} />
          </PieChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
